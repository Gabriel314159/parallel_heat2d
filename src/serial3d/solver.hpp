// The serial 3D pigment solver (model A of docs/design_3d.md): the 3D
// pigment in the 2D wind-driven flow of heat_serial.

#pragma once

#include <vector>

#include "common/io.hpp"
#include "common3d/config3d.hpp"

inline constexpr const char* serial3d_solver_name = "heat3d_serial";

struct Serial3DResult {
    RunSummary summary;
    std::vector<double> field;  // final field, nz x ny x nx, layer by layer (surface first), row by row
};

// Run the configured problem on one core. Writes snapshots if
// cfg.base.snapshot_interval > 0; the final field is returned, not written.
Serial3DResult run_serial_3d(const Config3D& c3);
