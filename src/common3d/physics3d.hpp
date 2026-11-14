// Initial conditions, the no-flux walls, mass diagnostics and the analytic
// cosine_decay reference of the 3D pigment (docs/design_3d.md, "Model A").
// The 3D counterparts of common/physics.hpp.

#pragma once

#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common3d/config3d.hpp"
#include "common3d/field3d.hpp"

// Depth coordinate of the centre of layer k; layers count down from the
// surface, so layer 1 is centred at z = lz - dz/2.
inline double cell_centre_z(int k, double dz, double lz) { return lz - (k - 0.5) * dz; }

// The initial value of the configured case at the point (x, y, z):
//   pigment_drop  c0 + A0 exp(-|r - r0|^2 / (2 sigma^2)), with r0 = (blob_x, blob_y, blob_depth);
//                 with --blob-z column, the 2D drop c0 + A0 exp(-((x-x0)^2 + (y-y0)^2) / (2 sigma^2))
//   constant      c0
//   cosine_decay  A0 cos(pi x/lx) cos(pi y/ly) cos(pi z/lz)
double initial_value_3d(const Config3D& c, double x, double y, double z);

// Set every interior cell of a local block (the x-y block of a rank, all
// layers) from its global coordinates, so every decomposition starts alike.
void set_initial_condition_3d(Field3D& field, const Config3D& c, const Block& b);

// No-flux walls: mirror the interior into the ghost cells across every wall
// of the block, the four side walls the block touches (as fill_wall_ghosts in
// 2D) and the surface and the bottom, which every block touches.
void fill_wall_ghosts_3d(Field3D& field, const Block& b);

// Interior sums of c dx dy dz and |c| dx dy dz for one block.
double local_mass_3d(const Field3D& field, double dx, double dy, double dz);
double local_abs_mass_3d(const Field3D& field, double dx, double dy, double dz);

// cosine_decay in 3D. Without flow, the mode is an exact eigenvector of the
// scheme: after n steps the field is (1 - lambda dt)^n times the initial one, with
//   lambda = 4 kappa (sin^2(pi dx/(2 lx))/dx^2 + sin^2(pi dy/(2 ly))/dy^2 + sin^2(pi dz/(2 lz))/dz^2).
double cosine_decay_continuous_3d(const Config3D& c, double x, double y, double z, double t);
double cosine_decay_discrete_factor_3d(const Config3D& c, double dt, int steps);

// Largest |c - reference| over one block's interior, divided by |A0|.
CosineErrors cosine_decay_errors_3d(const Field3D& field, const Config3D& c, const Block& b, double dt, int steps);
