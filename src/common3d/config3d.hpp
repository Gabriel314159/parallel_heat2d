// Run configuration of the 3D pigment model (heat3d_serial, heat3d_parallel):
// the 2D configuration plus the depth. See docs/design_3d.md, "Model A".
//
// Model A puts the pigment in three dimensions and keeps the 2D wind-driven
// flow, the same in every layer and with no vertical velocity. So every 2D
// option keeps its meaning (the x-y grid, the flow, the pigment parameters,
// output, threads), and four options are added: --nz, --lz, --blob-z and
// --freeze-flow. --decomp takes slabs or pencils.

#pragma once

#include <optional>
#include <string>

#include "common/config.hpp"

struct Config3D {
    Config base;  // all 2D options: the x-y grid, the flow, the pigment, output and parallel settings

    // Layers and depth. With the 2D defaults (256 x 256 cells on 1 x 1),
    // dz = 0.125 / 32 = 1/256 = dx: cubic cells.
    int nz = 32;
    double lz = 0.125;

    // Depth of the centre of the pigment_drop blob; empty means lz / 2.
    std::optional<double> blob_z;
    // --blob-z column: the 2D drop in every layer, so the field does not
    // change with depth (each layer then evolves exactly as the 2D model).
    bool blob_column = false;

    // Advance the flow only during the spin-up; afterwards the pigment moves
    // in that fixed flow and no Poisson solve is done. This makes the 3D
    // pigment the work that is timed (see docs/design_3d.md, "Model A").
    bool freeze_flow = false;
};

// The command-line help text.
std::string usage_3d(const std::string& program);

// Read the 3D options, translate --decomp slabs|pencils, and parse the rest
// with the 2D parse_args(). Throws std::invalid_argument as parse_args does,
// and for --decomp blocks, which needs the 3D flow (model B).
Config3D parse_args_3d(int argc, char** argv);

// validate() of the 2D options, and the 3D ones.
void validate_3d(const Config3D& c);

// "slabs" (y split: the 2D rows) or "pencils" (x and y split: the 2D blocks).
std::string decomp_name_3d(Decomp d);

// The depth of the blob centre: blob_z, or lz / 2.
double blob_depth(const Config3D& c);

// Cell sizes.
inline double dz_of(const Config3D& c) { return c.lz / c.nz; }

// The stability limits of docs/design_3d.md, "Model A", "Time step":
//   a = u_max/dx + v_max/dy           (no vertical velocity)
//   d = 2 kappa (1/dx^2 + 1/dy^2 + 1/dz^2)
//   n = nu (1/dx^2 + 1/dy^2)          (the flow is 2D)
// with dt (a + d) <= 1, dt a <= 0.5 and dt n <= 0.1. Only the pigment limit
// differs from 2D: the vertical diffusion makes it tighter.
double stable_dt_3d(const Config3D& c, double u_max, double v_max);

// Throws std::runtime_error naming the violated limit.
void check_time_step_3d(const Config3D& c, double dt, double u_max, double v_max);

// cfl times the largest stable step for speed max_speed, or the explicit dt
// if given and stable for that speed.
double time_step_3d(const Config3D& c);
