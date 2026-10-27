// Unit test for src/common/io.cpp: field files and snapshot coarsening.
// Built and run by CTest: ctest --test-dir build

#include "common/io.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

std::string read_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

}  // namespace

int main() {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "heat2d_test_io";
    std::filesystem::remove_all(dir);

    // A 4 x 6 field with values that are not exact in binary.
    const int ny = 4;
    const int nx = 6;
    std::vector<double> values;
    for (int k = 0; k < ny * nx; ++k) values.push_back(0.1 * k + 1.0 / 3.0);

    std::cout << "CSV\n";
    {
        const std::string path = (dir / "field.csv").string();
        write_field(path, values, ny, nx);
        std::ifstream in(path);
        std::vector<double> back;
        std::string line;
        int rows = 0;
        while (std::getline(in, line)) {
            ++rows;
            std::stringstream cells(line);
            std::string cell;
            while (std::getline(cells, cell, ',')) back.push_back(std::strtod(cell.c_str(), nullptr));
        }
        check(rows == ny, "one line per row");
        check(back == values, "17 significant digits read back exactly");
    }

    std::cout << "NumPy .npy\n";
    {
        const std::filesystem::path path = dir / "field.npy";
        write_field(path.string(), values, ny, nx);
        const std::string bytes = read_file(path);
        const std::size_t header = static_cast<unsigned char>(bytes[8]) + 256u * static_cast<unsigned char>(bytes[9]);
        check(bytes.compare(0, 6, "\x93NUMPY") == 0, "magic string");
        check((10 + header) % 64 == 0, "header padded to a multiple of 64 bytes");
        check(bytes.find("'shape': (4, 6)") != std::string::npos, "shape (ny, nx) in the header");
        check(bytes.size() == 10 + header + values.size() * sizeof(double), "raw doubles after the header");
    }

    std::cout << "Coarsening\n";
    {
        const std::vector<double> coarse = coarsen_field(values, ny, nx, 2);
        check(coarse.size() == 6, "4 x 6 with k = 2 gives 2 x 3");
        const double expected = (values[0] + values[1] + values[6] + values[7]) / 4;
        check(std::abs(coarse[0] - expected) < 1e-15, "first block is the average of its 4 cells");
        double fine_sum = 0.0;
        double coarse_sum = 0.0;
        for (double v : values) fine_sum += v;
        for (double v : coarse) coarse_sum += v;
        check(std::abs(coarse_sum * 4 - fine_sum) < 1e-12 * fine_sum, "mass kept: 4 x sum(coarse) = sum(fine)");
        bool rejected = false;
        try {
            coarsen_field(values, ny, nx, 4);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "k that does not divide nx is rejected");
    }

    std::filesystem::remove_all(dir);
    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
