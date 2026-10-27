// Function bodies for common/io.hpp.

#include "common/io.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifndef HEAT2D_BUILD_TYPE
#define HEAT2D_BUILD_TYPE "unknown"
#endif

namespace {

// Builds one JSON object, key by key. Numbers keep 17 significant digits;
// NaN and infinity become null, which JSON readers accept.
class JsonObject {
public:
    void add(const std::string& key, const std::string& value) { field(key) << '"' << escape(value) << '"'; }
    void add(const std::string& key, const char* value) { add(key, std::string(value)); }
    void add(const std::string& key, int value) { field(key) << value; }
    void add(const std::string& key, double value) {
        if (std::isfinite(value)) {
            field(key) << std::setprecision(17) << value;
        } else {
            field(key) << "null";
        }
    }
    std::string str() const { return "{" + body_.str() + "}"; }

private:
    std::ostringstream& field(const std::string& key) {
        if (!first_) body_ << ", ";
        first_ = false;
        body_ << '"' << escape(key) << "\": ";
        return body_;
    }
    static std::string escape(const std::string& text) {
        std::string out;
        for (char ch : text) {
            if (ch == '"' || ch == '\\') {
                out += '\\';
                out += ch;
            } else if (ch == '\n') {
                out += "\\n";
            } else {
                out += ch;
            }
        }
        return out;
    }

    std::ostringstream body_;
    bool first_ = true;
};

void add_config(JsonObject& json, const Config& c) {
    json.add("nx", c.nx);
    json.add("ny", c.ny);
    json.add("lx", c.lx);
    json.add("ly", c.ly);
    json.add("steps", c.steps);
    json.add("kappa", c.kappa);
    json.add("cfl", c.cfl);
    json.add("flow_amplitude", c.flow_amplitude);
    json.add("flow_epsilon", c.flow_epsilon);
    json.add("flow_period", c.flow_period);
    json.add("case", case_name(c.initial_case));
    json.add("background", c.background);
    json.add("amplitude", c.amplitude);
    json.add("blob_x", c.blob_x);
    json.add("blob_y", c.blob_y);
    json.add("blob_width", c.blob_width);
}

std::ofstream open_for_writing(const std::string& path, std::ios::openmode mode = std::ios::out) {
    ensure_parent_directory(path);
    std::ofstream out(path, mode);
    if (!out) throw std::runtime_error("cannot open " + path + " for writing");
    return out;
}

}  // namespace

std::string compiler_id() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("gcc ") + __VERSION__;
#else
    return "unknown";
#endif
}

std::string build_type() {
    return HEAT2D_BUILD_TYPE;
}

std::string summary_json(const Config& c, const RunSummary& s) {
    JsonObject json;
    json.add("solver", s.solver);
    add_config(json, c);
    json.add("dt", s.dt);
    json.add("t_end", s.t_end);
    json.add("mpi_ranks", s.mpi_ranks);
    json.add("omp_threads", s.omp_threads);
    json.add("decomposition", s.decomposition);
    json.add("exchange", s.exchange);
    json.add("py", s.py);
    json.add("px", s.px);
    json.add("elapsed_s", s.elapsed_s);
    json.add("snapshot_s", s.snapshot_s);
    const double updates = static_cast<double>(c.nx) * c.ny * c.steps;
    json.add("cell_updates_per_s", s.elapsed_s > 0 ? updates / s.elapsed_s : 0.0);
    json.add("mass_initial", s.mass_initial);
    json.add("mass_final", s.mass_final);
    json.add("mass_error", s.mass_error);
    json.add("peclet", s.peclet);
    json.add("peclet_grid", s.peclet_grid);
    json.add("cosine_discrete_error", s.cosine_discrete_error);
    json.add("cosine_continuous_error", s.cosine_continuous_error);
    json.add("compiler", s.compiler);
    json.add("mpi_library", s.mpi_library);
    json.add("build_type", s.build_type);
    return json.str();
}

void print_summary(std::ostream& out, const Config& c, const RunSummary& s) {
    const double updates = static_cast<double>(c.nx) * c.ny * c.steps;
    out << s.solver << ": " << c.nx << " x " << c.ny << " cells, " << c.steps << " steps, case "
        << case_name(c.initial_case) << "\n"
        << "  dt = " << s.dt << ", final time = " << s.t_end << "\n"
        << "  flow A = " << c.flow_amplitude << ", epsilon = " << c.flow_epsilon << ", period = " << c.flow_period
        << ", kappa = " << c.kappa << "\n"
        << "  Peclet = " << s.peclet << ", grid Peclet = " << s.peclet_grid << "\n"
        << "  solver loop: " << s.elapsed_s << " s ("
        << (s.elapsed_s > 0 ? updates / s.elapsed_s : 0.0) << " cell updates/s)";
    if (s.snapshot_s > 0) out << ", snapshots: " << s.snapshot_s << " s";
    out << "\n  mass: initial " << s.mass_initial << ", final " << s.mass_final
        << ", relative change " << s.mass_error << "\n";
    if (std::isfinite(s.cosine_discrete_error)) {
        out << "  cosine_decay error: vs exact discrete " << s.cosine_discrete_error << ", vs continuous "
            << s.cosine_continuous_error << "\n";
    }
    out << "RESULT " << summary_json(c, s) << "\n";
}

std::vector<double> copy_interior(const Field& field) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(field.ny()) * static_cast<std::size_t>(field.nx()));
    for (int i = 1; i <= field.ny(); ++i) {
        const double* row = field.row(i);
        values.insert(values.end(), row + 1, row + 1 + field.nx());
    }
    return values;
}

void ensure_parent_directory(const std::string& path) {
    const std::filesystem::path parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}

std::string replace_extension(const std::string& path, const std::string& extension) {
    return std::filesystem::path(path).replace_extension(extension).string();
}

void write_field_csv(const std::string& path, const std::vector<double>& values, int ny, int nx) {
    if (values.size() != static_cast<std::size_t>(ny) * static_cast<std::size_t>(nx)) {
        throw std::invalid_argument("write_field_csv: size does not match ny x nx");
    }
    std::ofstream out = open_for_writing(path);
    char number[32];
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            std::snprintf(number, sizeof number, "%.17g", values[static_cast<std::size_t>(i) * nx + j]);
            if (j > 0) out << ',';
            out << number;
        }
        out << '\n';
    }
    if (!out) throw std::runtime_error("error while writing " + path);
}

// NumPy .npy version 1.0: magic, version, header length, a Python-dict header
// padded to a multiple of 64 bytes, then the raw little-endian doubles.
void write_field_npy(const std::string& path, const std::vector<double>& values, int ny, int nx) {
    if (values.size() != static_cast<std::size_t>(ny) * static_cast<std::size_t>(nx)) {
        throw std::invalid_argument("write_field_npy: size does not match ny x nx");
    }
    const std::uint16_t probe = 1;
    if (*reinterpret_cast<const unsigned char*>(&probe) != 1) {
        throw std::runtime_error("write_field_npy assumes a little-endian machine");
    }

    std::string header = "{'descr': '<f8', 'fortran_order': False, 'shape': (" + std::to_string(ny) + ", " +
                         std::to_string(nx) + "), }";
    const std::size_t preamble = 10;  // magic (6) + version (2) + header length (2)
    const std::size_t unpadded = preamble + header.size() + 1;
    header.append((unpadded + 63) / 64 * 64 - unpadded, ' ');
    header.push_back('\n');

    std::ofstream out = open_for_writing(path, std::ios::out | std::ios::binary);
    out.write("\x93NUMPY", 6);
    const char version[2] = {1, 0};
    out.write(version, 2);
    const char length[2] = {static_cast<char>(header.size() & 0xff), static_cast<char>((header.size() >> 8) & 0xff)};
    out.write(length, 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(double)));
    if (!out) throw std::runtime_error("error while writing " + path);
}

void write_text_file(const std::string& path, const std::string& text) {
    std::ofstream out = open_for_writing(path);
    out << text;
    if (!out) throw std::runtime_error("error while writing " + path);
}

void write_field(const std::string& path, const std::vector<double>& values, int ny, int nx) {
    if (std::filesystem::path(path).extension() == ".npy") {
        write_field_npy(path, values, ny, nx);
    } else {
        write_field_csv(path, values, ny, nx);
    }
}

std::vector<double> coarsen_field(const std::vector<double>& values, int ny, int nx, int k) {
    if (k < 1 || ny % k != 0 || nx % k != 0) {
        throw std::invalid_argument("coarsen_field: k must divide ny and nx");
    }
    const int cy = ny / k;
    const int cx = nx / k;
    std::vector<double> coarse(static_cast<std::size_t>(cy) * static_cast<std::size_t>(cx), 0.0);
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            coarse[static_cast<std::size_t>(i / k) * cx + j / k] += values[static_cast<std::size_t>(i) * nx + j];
        }
    }
    const double cells = static_cast<double>(k) * k;
    for (double& value : coarse) value /= cells;
    return coarse;
}

SnapshotWriter::SnapshotWriter(const Config& c, double dt, const std::string& solver)
    : dir_(c.snapshot_dir), dt_(dt), ny_(c.ny), nx_(c.nx), coarsen_(c.snapshot_coarsen) {
    std::filesystem::create_directories(dir_);

    JsonObject meta;
    meta.add("solver", solver);
    add_config(meta, c);
    meta.add("dt", dt);
    meta.add("snapshot_interval", c.snapshot_interval);
    meta.add("snapshot_coarsen", coarsen_);
    meta.add("snapshot_ny", ny_ / coarsen_);
    meta.add("snapshot_nx", nx_ / coarsen_);
    write_text_file(dir_ + "/meta.json", meta.str() + "\n");

    index_ = open_for_writing(dir_ + "/index.csv");
    index_ << "step,time,mass,file\n";
}

void SnapshotWriter::write(int step, const std::vector<double>& values, double mass) {
    char name[32];
    std::snprintf(name, sizeof name, "field_%06d.npy", step);
    if (coarsen_ > 1) {
        write_field_npy(dir_ + "/" + name, coarsen_field(values, ny_, nx_, coarsen_), ny_ / coarsen_, nx_ / coarsen_);
    } else {
        write_field_npy(dir_ + "/" + name, values, ny_, nx_);
    }
    index_ << step << ',' << std::setprecision(17) << step * dt_ << ',' << mass << ',' << name << '\n';
    index_.flush();
    if (!index_) throw std::runtime_error("error while writing " + dir_ + "/index.csv");
}
