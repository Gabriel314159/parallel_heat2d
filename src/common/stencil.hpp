// The conservative upwind / central-difference update of one cell.
// See docs/design.md, "Discretization", and src/README.md, "Stencil".
//
// No communication, threading, allocation or output happens here: the serial
// and parallel solvers decide which rows to update and call update_row().

#pragma once

#include "common/grid.hpp"

// Constants of the update, computed once per run.
struct StencilCoeffs {
    double kappa_dx;  // kappa / dx
    double kappa_dy;  // kappa / dy
    double dt_dx;     // dt / dx
    double dt_dy;     // dt / dy
};

inline StencilCoeffs make_stencil_coeffs(double kappa, double dx, double dy, double dt) {
    return StencilCoeffs{kappa / dx, kappa / dy, dt / dx, dt / dy};
}

// Total (advective + diffusive) flux through one face, counted positive from
// the "minus" cell to the "plus" cell: west -> east in x, south -> north in y.
// vel is the face velocity in the same direction, kappa_h is kappa / dx or
// kappa / dy. The upwind value is taken from the cell the flow comes from.
//
// Both cells sharing a face call this with the same arguments in the same
// order, so they get the identical flux and the total mass is conserved.
inline double face_flux(double c_minus, double c_plus, double vel, double kappa_h) {
    const double upwind = vel > 0.0 ? c_minus : c_plus;
    return vel * upwind - kappa_h * (c_plus - c_minus);
}

// New value of one cell from its own value, its four neighbours and the
// velocities on its four faces (u: x-velocity, positive east; v: y-velocity,
// positive north).
inline double update_cell(double c, double c_west, double c_east, double c_north, double c_south,
                          double u_west, double u_east, double v_north, double v_south,
                          const StencilCoeffs& k) {
    const double f_east = face_flux(c, c_east, u_east, k.kappa_dx);
    const double f_west = face_flux(c_west, c, u_west, k.kappa_dx);
    const double g_north = face_flux(c, c_north, v_north, k.kappa_dy);
    const double g_south = face_flux(c_south, c, v_south, k.kappa_dy);
    return c - k.dt_dx * (f_east - f_west) - k.dt_dy * (g_north - g_south);
}

// Update cells (i, j_first..j_last) of `next` from `c`. Rows are numbered
// top-down, so the north neighbour of row i is row i-1.
inline void update_row(const Field& c, Field& next, const FaceVelocities& vel, const StencilCoeffs& k,
                       int i, int j_first, int j_last) {
    const double* mid = c.row(i);
    const double* above = c.row(i - 1);
    const double* below = c.row(i + 1);
    const double* u = vel.east_row(i);         // u[j] = east(i, j); west face of (i, j) is u[j - 1]
    const double* v_n = vel.north_row(i);      // v_n[j - 1] = north(i, j)
    const double* v_s = vel.north_row(i + 1);  // v_s[j - 1] = north(i + 1, j) = south face of (i, j)
    double* out = next.row(i);
    for (int j = j_first; j <= j_last; ++j) {
        out[j] = update_cell(mid[j], mid[j - 1], mid[j + 1], above[j], below[j],
                             u[j - 1], u[j], v_n[j - 1], v_s[j - 1], k);
    }
}
