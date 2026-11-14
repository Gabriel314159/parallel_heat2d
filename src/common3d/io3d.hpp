// Output of the 3D pigment model (docs/design_3d.md, "Model A", "Output").
//
// - Final field (--output PATH): the whole 3D field, as CSV with nz x ny rows
//   of nx values (layer by layer, the surface layer first; within a layer,
//   the north row first, as in 2D) or as .npy of shape (nz, ny, nx). The run
//   summary goes next to it as .json, as in 2D.
// - Snapshots: the 2D snapshot folder of common/io.hpp, whose pigment is the
//   surface layer (so the 2D plotting scripts work unchanged), plus
//     section_NNNNNN.npy  a vertical section, nz x nx: the cells of global
//                         row section_row in every layer, surface first,
//     meta3d.json         nz, lz, dz, blob depth and the section's row and y.
// - Run summary: the 2D RESULT line with the 3D fields added (model, nz, lz,
//   blob_z, flow, pigment_s, pigment_halo_s, section_row).

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "common/io.hpp"
#include "common3d/config3d.hpp"
#include "common3d/field3d.hpp"

// Interior values of a block, nz x ny x nx, layer by layer, row by row.
std::vector<double> copy_interior_3d(const Field3D& field);

// Interior of one layer k, ny x nx.
std::vector<double> copy_layer(const Field3D& field, int k);

// Row i of every layer, nz x nx (the interior columns).
std::vector<double> copy_section(const Field3D& field, int i);

// CSV (nz ny rows of nx values) or .npy of shape (nz, ny, nx), by extension.
void write_field_3d(const std::string& path, const std::vector<double>& values, int nz, int ny, int nx);

// The global row of the vertical section: the row whose cells contain y = blob_y.
int section_row(const Config3D& c);

// What the 3D runs report in addition to the 2D summary.
struct Report3D {
    double pigment_s = 0.0;       // pigment updates of the main loop (wall ghosts and halos included); max over ranks
    double pigment_halo_s = 0.0;  // their halo exchanges (posting and waiting, when overlapping); max over ranks
};

// Fill the 3D parts of a summary: the grid text, cells per step, the extra
// RESULT fields and the report lines. Call it after the 2D fields are set.
void describe_3d(const Config3D& c, const Report3D& r, RunSummary& s);

// The snapshot folder: SnapshotWriter for the surface layer, plus sections.
class SnapshotWriter3D {
public:
    SnapshotWriter3D(const Config3D& c, double dt, double t_release, const std::string& solver);
    // `surface` holds the surface layer as its pigment; `section` is nz x nx.
    void write(const Snapshot& surface, const std::vector<double>& section);

private:
    SnapshotWriter surface_;
    std::string dir_;
    int nz_;
    int nx_;
};
