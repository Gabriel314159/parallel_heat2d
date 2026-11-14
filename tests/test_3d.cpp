// Unit tests for the serial 3D pigment model (src/common3d, src/serial3d):
// model A of docs/design_3d.md, the 3D pigment in the 2D flow.
// Built and run by CTest: ctest --test-dir build
//
// The strongest checks compare with the 2D solver, which is already
// validated: with one layer, or with a field that does not change with depth,
// every layer must evolve exactly as the 2D pigment.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "common3d/config3d.hpp"
#include "common3d/field3d.hpp"
#include "common3d/io3d.hpp"
#include "common3d/physics3d.hpp"
#include "common3d/stencil3d.hpp"
#include "serial/solver.hpp"
#include "serial3d/solver.hpp"

namespace {

int failures = 0;

std::string sci(double value) {
    std::ostringstream out;
    out.precision(2);
    out << std::scientific << value;
    return out.str();
}

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

// The small basin of tests/test_serial_solver.cpp: 1 x 1, a developed flow at
// Re = 50 after a short spin-up, and a blob several cells wide.
Config small_basin(int nx, int ny, int steps) {
    Config c;
    c.nx = nx;
    c.ny = ny;
    c.steps = steps;
    c.spinup = 1.0;
    c.viscosity = 0.02;
    c.wind_amplitude = 0.02 * std::pow(2 * M_PI / c.ly, 3);
    c.max_speed = 2.0;
    c.blob_width = 0.15;
    c.blob_y = 0.5;
    c.blob_x = 0.3;
    c.check_interval = 10;
    return c;
}

// Parse a command line as heat3d_serial would.
Config3D parse(std::vector<std::string> words) {
    words.insert(words.begin(), "heat3d_serial");
    std::vector<char*> argv;
    for (std::string& w : words) argv.push_back(w.data());
    return parse_args_3d(static_cast<int>(argv.size()), argv.data());
}

// Largest |a - b| over two fields of the same size.
double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return INFINITY;
    double worst = 0.0;
    for (std::size_t n = 0; n < a.size(); ++n) worst = std::max(worst, std::abs(a[n] - b[n]));
    return worst;
}

}  // namespace

int main() {
    std::cout << "Field3D storage\n";
    {
        Field3D f(3, 4, 5);
        check(f.size() == static_cast<std::size_t>(5 * 6 * 7), "(nz+2)(ny+2)(nx+2) = 210 doubles, ghosts included");
        check(f.layer_stride() == static_cast<std::size_t>(6 * 7) && f.stride() == 7, "row and layer strides");
        bool zero = true;
        for (double v : f.values()) zero = zero && v == 0.0;
        check(zero, "starts at 0.0");
        f(2, 3, 4) = 7.5;
        check(f.row(2, 3)[4] == 7.5 && &f(2, 3, 5) == &f(2, 3, 4) + 1, "row(k, i)[j] is (k, i, j); x is contiguous");
        check(&f(3, 3, 4) == &f(2, 3, 4) + f.layer_stride(), "one layer down is one layer stride further");
    }

    std::cout << "Command line\n";
    {
        const Config3D d = parse({});
        check(d.nz == 32 && d.lz == 0.125 && !d.freeze_flow && !d.blob_column && !d.blob_z,
              "defaults: 32 layers, depth 0.125, live flow, blob at lz/2");
        check(std::abs(dz_of(d) - d.base.lx / d.base.nx) < 1e-15, "default cells are cubes (dz = dx)");
        const Config3D c = parse({"--nz", "8", "--lz", "0.5", "--blob-z", "0.1", "--freeze-flow", "--decomp", "pencils",
                                  "--nx", "40", "--kappa", "0.01"});
        validate_3d(c);
        check(c.nz == 8 && c.lz == 0.5 && c.blob_z && *c.blob_z == 0.1 && c.freeze_flow, "--nz, --lz, --blob-z, --freeze-flow");
        check(c.base.decomposition == Decomp::two_d && decomp_name_3d(c.base.decomposition) == "pencils",
              "--decomp pencils is the 2D blocks");
        check(c.base.nx == 40 && c.base.kappa == 0.01, "the 2D options pass through");
        check(parse({"--decomp", "slabs"}).base.decomposition == Decomp::rows, "--decomp slabs is the 2D rows");
        check(parse({"--blob-z", "column"}).blob_column, "--blob-z column");
        bool blocks = false;
        try {
            parse({"--decomp", "blocks"});
        } catch (const std::invalid_argument&) {
            blocks = true;
        }
        check(blocks, "--decomp blocks is refused (it needs the 3D flow)");
        bool bad_nz = false;
        try {
            validate_3d(parse({"--nz", "0"}));
        } catch (const std::invalid_argument&) {
            bad_nz = true;
        }
        check(bad_nz, "--nz 0 is refused");
    }

    std::cout << "Time step\n";
    {
        Config3D c;
        c.base = small_basin(48, 48, 1);
        c.nz = 4;
        c.lz = 0.25;
        const double dt3 = time_step_3d(c);
        check(dt3 <= time_step(c.base) * (1 + 1e-15), "never larger than the 2D step");
        // Thin layers: vertical diffusion decides the pigment limit.
        c.nz = 400;
        c.base.kappa = 0.01;
        const double dx = 1.0 / 48, dz = dz_of(c);
        const double a = c.base.max_speed / dx + c.base.max_speed / dx;
        const double d = 2 * c.base.kappa * (2 / (dx * dx) + 1 / (dz * dz));
        check(std::abs(stable_dt_3d(c, c.base.max_speed, c.base.max_speed) - 1 / (a + d)) < 1e-15 / (a + d),
              "with thin layers the limit is 1 / (a + d), d with the 1/dz^2 term");
        bool rejected = false;
        try {
            check_time_step_3d(c, 2 * stable_dt_3d(c, 1.0, 1.0), 1.0, 1.0);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check(rejected, "check_time_step_3d rejects twice the stable dt");
    }

    std::cout << "One layer is the 2D model\n";
    {
        const Config base = small_basin(37, 29, 60);
        const SerialResult r2 = run_serial(base);
        Config3D c;
        c.base = base;
        c.base.dt = r2.summary.dt;  // the same step: the 3D limit alone would choose a smaller one
        c.nz = 1;
        c.lz = 0.4;
        const Serial3DResult r3 = run_serial_3d(c);
        check(r3.field.size() == r2.field.size() && max_difference(r3.field, r2.field) == 0.0,
              "the pigment equals the 2D pigment exactly (max difference " + sci(max_difference(r3.field, r2.field)) + ")");
        check(r3.summary.energy == r2.summary.energy && r3.summary.u_max_measured == r2.summary.u_max_measured,
              "the flow is the 2D flow exactly");
        check(r3.summary.mass_error < 1e-12, "mass conserved (" + sci(r3.summary.mass_error) + ")");
    }

    std::cout << "A field that does not change with depth stays so: every layer is the 2D model\n";
    {
        const Config base = small_basin(37, 29, 60);
        const SerialResult r2 = run_serial(base);
        Config3D c;
        c.base = base;
        c.base.dt = r2.summary.dt;
        c.nz = 5;
        c.lz = 0.3;
        c.blob_column = true;
        const Serial3DResult r3 = run_serial_3d(c);
        const std::size_t layer = r2.field.size();
        double worst = 0.0;
        for (int k = 0; k < c.nz; ++k) {
            const std::vector<double> one(r3.field.begin() + static_cast<std::ptrdiff_t>(k * layer),
                                          r3.field.begin() + static_cast<std::ptrdiff_t>((k + 1) * layer));
            worst = std::max(worst, max_difference(one, r2.field));
        }
        check(worst == 0.0, "each of the 5 layers equals the 2D pigment exactly (max difference " + sci(worst) + ")");
    }

    std::cout << "3D blob in the flow\n";
    {
        Config3D c;
        c.base = small_basin(32, 32, 120);
        c.nz = 12;
        c.lz = 0.375;  // dz = dx
        c.base.blob_width = 0.08;
        const Serial3DResult r = run_serial_3d(c);
        check(r.summary.mass_error < 1e-12, "relative mass change < 1e-12 (" + sci(r.summary.mass_error) + ")");
        const double lowest = *std::min_element(r.field.begin(), r.field.end());
        const double highest = *std::max_element(r.field.begin(), r.field.end());
        check(lowest >= -1e-14 && highest <= 1.0 + 1e-14, "no new extrema: the pigment stays within [0, 1]");
        // A blob centred in depth, in a flow without vertical velocity and
        // between two mirrored walls, stays symmetric about mid-depth.
        const std::size_t layer = static_cast<std::size_t>(c.base.nx) * c.base.ny;
        double asymmetry = 0.0;
        for (int k = 0; k < c.nz / 2; ++k) {
            for (std::size_t n = 0; n < layer; ++n) {
                asymmetry = std::max(asymmetry, std::abs(r.field[k * layer + n] - r.field[(c.nz - 1 - k) * layer + n]));
            }
        }
        check(asymmetry < 1e-14, "symmetric about mid-depth (" + sci(asymmetry) + ")");
        check(r.summary.cells == 32.0 * 32 * 12 && r.summary.extra_json.find("\"nz\": 12") != std::string::npos,
              "the summary reports the 3D grid");
    }

    std::cout << "Constant field under the flow\n";
    {
        Config3D c;
        c.base = small_basin(24, 24, 80);
        c.base.initial_case = InitialCase::constant;
        c.base.background = 42.0;
        c.nz = 6;
        c.lz = 0.25;
        const Serial3DResult r = run_serial_3d(c);
        double deviation = 0.0;
        for (double v : r.field) deviation = std::max(deviation, std::abs(v - 42.0));
        check(deviation <= 1e-12 * 42.0, "stays 42 to rounding (max deviation " + sci(deviation) + ")");
    }

    std::cout << "Frozen flow\n";
    {
        Config3D c;
        c.base = small_basin(24, 24, 50);
        c.nz = 4;
        c.lz = 0.2;
        c.freeze_flow = true;
        const Serial3DResult frozen = run_serial_3d(c);
        check(frozen.summary.poisson_solves == 0 && frozen.summary.poisson_iterations == 0,
              "no Poisson solve after the spin-up");
        check(frozen.summary.u_max_measured > 0.1, "the pigment moves in the spun-up flow");
        check(frozen.summary.mass_error < 1e-12, "mass conserved (" + sci(frozen.summary.mass_error) + ")");
        check(frozen.summary.extra_json.find("\"flow\": \"frozen\"") != std::string::npos, "reported as frozen");
    }

    std::cout << "Cosine decay without wind\n";
    {
        Config3D c;
        c.base.nx = 20;
        c.base.ny = 16;
        c.base.steps = 40;
        c.base.initial_case = InitialCase::cosine_decay;
        c.base.wind_amplitude = 0.0;
        c.base.max_speed = 0.0;
        c.base.spinup = 0.0;
        c.base.kappa = 1e-3;
        c.nz = 12;
        c.lz = 0.5;
        const Serial3DResult r = run_serial_3d(c);
        check(r.summary.cosine_discrete_error < 1e-12,
              "matches (1 - lambda dt)^n exactly (" + sci(r.summary.cosine_discrete_error) + ")");
        check(r.summary.mass_error < 1e-12, "zero total mass stays zero");
    }

    std::cout << "Cosine decay converges under grid refinement\n";
    {
        double previous = 0.0;
        bool decreasing = true;
        for (int n : {8, 16, 32}) {
            Config3D c;
            c.base.nx = n;
            c.base.ny = n;
            c.base.steps = 1;
            c.base.initial_case = InitialCase::cosine_decay;
            c.base.wind_amplitude = 0.0;
            c.base.max_speed = 0.0;
            c.base.spinup = 0.0;
            c.base.kappa = 0.005;
            c.base.viscosity = 1e-9;  // keeps the flow's limit out of the way
            c.base.cfl = 0.5;
            c.nz = n / 2;
            c.lz = 0.5;
            const double dt = time_step_3d(c);
            c.base.steps = static_cast<int>(std::lround(0.5 / dt));
            const Serial3DResult r = run_serial_3d(c);
            std::cout << "    " << n << " x " << n << " x " << n / 2
                      << ": error vs continuous = " << sci(r.summary.cosine_continuous_error) << "\n";
            if (previous > 0) decreasing = decreasing && r.summary.cosine_continuous_error < previous / 3;
            previous = r.summary.cosine_continuous_error;
        }
        check(decreasing, "error drops more than 3x per refinement (at least second order)");
    }

    std::cout << "Output\n";
    {
        Field3D f(2, 3, 4);
        for (int k = 1; k <= 2; ++k)
            for (int i = 1; i <= 3; ++i)
                for (int j = 1; j <= 4; ++j) f(k, i, j) = 100 * k + 10 * i + j;
        const std::vector<double> all = copy_interior_3d(f);
        check(all.size() == 24 && all[0] == 111 && all[23] == 234, "whole field: layer by layer, row by row");
        const std::vector<double> top = copy_layer(f, 1);
        const std::vector<double> section = copy_section(f, 2);
        check(top.size() == 12 && top[11] == 134, "surface layer");
        check(section.size() == 8 && section[0] == 121 && section[7] == 224, "vertical section: row 2 of every layer");
        Config3D c;
        c.base.ny = 10;
        c.base.ly = 1.0;
        c.base.blob_y = 0.55;
        check(section_row(c) == 5, "the section row contains y = blob_y (row 5 spans 0.5 to 0.6)");
    }

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
