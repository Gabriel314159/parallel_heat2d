// Function bodies for common3d/io3d.hpp.

#include "common3d/io3d.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <stdexcept>

std::vector<double> copy_interior_3d(const Field3D& field) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(field.nz()) * static_cast<std::size_t>(field.ny()) *
                   static_cast<std::size_t>(field.nx()));
    for (int k = 1; k <= field.nz(); ++k) {
        for (int i = 1; i <= field.ny(); ++i) {
            const double* row = field.row(k, i);
            values.insert(values.end(), row + 1, row + 1 + field.nx());
        }
    }
    return values;
}

std::vector<double> copy_layer(const Field3D& field, int k) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(field.ny()) * static_cast<std::size_t>(field.nx()));
    for (int i = 1; i <= field.ny(); ++i) {
        const double* row = field.row(k, i);
        values.insert(values.end(), row + 1, row + 1 + field.nx());
    }
    return values;
}

std::vector<double> copy_section(const Field3D& field, int i) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(field.nz()) * static_cast<std::size_t>(field.nx()));
    for (int k = 1; k <= field.nz(); ++k) {
        const double* row = field.row(k, i);
        values.insert(values.end(), row + 1, row + 1 + field.nx());
    }
    return values;
}

void write_field_3d(const std::string& path, const std::vector<double>& values, int nz, int ny, int nx) {
    if (std::filesystem::path(path).extension() == ".npy") {
        write_npy(path, values,
                  {static_cast<std::size_t>(nz), static_cast<std::size_t>(ny), static_cast<std::size_t>(nx)});
    } else {
        write_field_csv(path, values, nz * ny, nx);
    }
}

int section_row(const Config3D& c) {
    const double dy = c.base.ly / c.base.ny;
    const int row = static_cast<int>((c.base.ly - c.base.blob_y) / dy) + 1;
    return std::clamp(row, 1, c.base.ny);
}

void describe_3d(const Config3D& c, const Report3D& r, RunSummary& s) {
    const Config& b = c.base;
    s.grid = std::to_string(b.nx) + " x " + std::to_string(b.ny) + " x " + std::to_string(c.nz);
    s.cells = static_cast<double>(b.nx) * b.ny * c.nz;
    s.decomposition = s.decomposition == "none" ? "none" : decomp_name_3d(b.decomposition);

    JsonObject extra;
    extra.add("model", "3d_pigment_2d_flow");
    extra.add("nz", c.nz);
    extra.add("lz", c.lz);
    extra.add("dz", dz_of(c));
    if (c.blob_column) {
        extra.add("blob_z", "column");
    } else {
        extra.add("blob_z", blob_depth(c));
    }
    extra.add("flow", c.freeze_flow ? "frozen" : "live");
    extra.add("pigment_s", r.pigment_s);
    extra.add("pigment_halo_s", r.pigment_halo_s);
    extra.add("section_row", section_row(c));
    s.extra_json = extra.fields();

    std::ostringstream lines;
    lines << "  3D pigment: " << c.nz << " layers over depth " << c.lz << " (dz = " << dz_of(c) << "), ";
    if (b.initial_case == InitialCase::pigment_drop) {
        if (c.blob_column) {
            lines << "the 2D drop in every layer, ";
        } else {
            lines << "blob centre at depth z = " << blob_depth(c) << ", ";
        }
    }
    lines << (c.freeze_flow ? "flow frozen after the spin-up" : "flow advanced every step") << "\n"
          << "  pigment updates: " << r.pigment_s << " s";
    if (s.solver.find("parallel") != std::string::npos) lines << " (halo exchanges " << r.pigment_halo_s << " s)";
    lines << "\n";
    s.report_lines = lines.str();
}

SnapshotWriter3D::SnapshotWriter3D(const Config3D& c, double dt, double t_release, const std::string& solver)
    : surface_(c.base, dt, t_release, solver), dir_(c.base.snapshot_dir), nz_(c.nz), nx_(c.base.nx) {
    const int row = section_row(c);
    JsonObject meta;
    meta.add("solver", solver);
    meta.add("model", "3d_pigment_2d_flow");
    meta.add("nz", c.nz);
    meta.add("lz", c.lz);
    meta.add("dz", dz_of(c));
    meta.add("lx", c.base.lx);
    meta.add("nx", c.base.nx);
    meta.add("blob_z", c.blob_column ? -1.0 : blob_depth(c));
    meta.add("blob_column", c.blob_column);
    meta.add("flow", c.freeze_flow ? "frozen" : "live");
    meta.add("section_row", row);
    meta.add("section_y", c.base.ly - (row - 0.5) * c.base.ly / c.base.ny);
    meta.add("surface_layer", "field_NNNNNN.npy is layer 1, the top layer");
    write_text_file(dir_ + "/meta3d.json", meta.str() + "\n");
}

void SnapshotWriter3D::write(const Snapshot& surface, const std::vector<double>& section) {
    surface_.write(surface);
    char name[40];
    std::snprintf(name, sizeof name, "section_%06d.npy", surface.step);
    write_npy(dir_ + "/" + name, section, {static_cast<std::size_t>(nz_), static_cast<std::size_t>(nx_)});
}
