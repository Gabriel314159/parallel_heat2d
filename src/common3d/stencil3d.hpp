// The conservative upwind / central-difference update of one 3D cell
// (docs/design_3d.md, "Model A", "Pigment update").
//
// The 2D update of common/stencil.hpp with a third pair of faces, the top and
// bottom of the cell. The flux through every face is face_flux() of the 2D
// code, positive along +x, +y and +z (up). Both cells sharing a face compute
// that flux with the same arguments in the same order, so the total pigment
// is conserved to rounding, as in 2D.
//
// Model A has no vertical velocity: the flow is the 2D flow, the same in
// every layer, so w = 0 on every top and bottom face and only diffusion
// crosses them. update_cell_3d() still takes w, so that a 3D flow can use it.
//
// No communication or allocation happens here; update_region_3d shares the
// (layer, row) pairs out among the OpenMP threads when the caller is built
// with OpenMP.

#pragma once

#include "common/grid.hpp"
#include "common/stencil.hpp"
#include "common3d/field3d.hpp"

struct StencilCoeffs3D {
    double kappa_dx;  // kappa / dx
    double kappa_dy;  // kappa / dy
    double kappa_dz;  // kappa / dz
    double dt_dx;     // dt / dx
    double dt_dy;     // dt / dy
    double dt_dz;     // dt / dz
};

inline StencilCoeffs3D make_stencil_coeffs_3d(double kappa, double dx, double dy, double dz, double dt) {
    return StencilCoeffs3D{kappa / dx, kappa / dy, kappa / dz, dt / dx, dt / dy, dt / dz};
}

// New value of one cell from its own value, its six neighbours and the
// velocities on its six faces. "up" is the layer above (k - 1, larger z),
// "down" the layer below (k + 1).
//
// With c_up == c_down == c (a field that does not change with depth, or a
// single layer between two mirrored wall ghosts), both vertical fluxes are
// exactly zero, and the result equals the 2D update_cell() exactly (only the
// sign of a result of zero could differ).
inline double update_cell_3d(double c, double c_west, double c_east, double c_north, double c_south, double c_up,
                             double c_down, double u_west, double u_east, double v_north, double v_south, double w_up,
                             double w_down, const StencilCoeffs3D& k) {
    const double f_east = face_flux(c, c_east, u_east, k.kappa_dx);
    const double f_west = face_flux(c_west, c, u_west, k.kappa_dx);
    const double g_north = face_flux(c, c_north, v_north, k.kappa_dy);
    const double g_south = face_flux(c_south, c, v_south, k.kappa_dy);
    const double h_top = face_flux(c, c_up, w_up, k.kappa_dz);         // minus: this cell, plus: the one above
    const double h_bottom = face_flux(c_down, c, w_down, k.kappa_dz);  // minus: the one below, plus: this cell
    return c - k.dt_dx * (f_east - f_west) - k.dt_dy * (g_north - g_south) - k.dt_dz * (h_top - h_bottom);
}

// Update cells (layer, i, j_first..j_last) of `next` from `c`, with the 2D
// face velocities of row i (the same in every layer) and w = 0.
inline void update_row_3d(const Field3D& c, Field3D& next, const FaceVelocities& vel, const StencilCoeffs3D& k,
                          int layer, int i, int j_first, int j_last) {
    const double* mid = c.row(layer, i);
    const double* above = c.row(layer, i - 1);  // north neighbour
    const double* below = c.row(layer, i + 1);  // south neighbour
    const double* up = c.row(layer - 1, i);     // the layer above
    const double* down = c.row(layer + 1, i);   // the layer below
    const double* u = vel.east_row(i);
    const double* v_n = vel.north_row(i);
    const double* v_s = vel.north_row(i + 1);
    double* out = next.row(layer, i);
    for (int j = j_first; j <= j_last; ++j) {
        out[j] = update_cell_3d(mid[j], mid[j - 1], mid[j + 1], above[j], below[j], up[j], down[j], u[j - 1], u[j],
                                v_n[j - 1], v_s[j - 1], 0.0, 0.0, k);
    }
}

// Update rows i0..i1, columns j0..j1 of every layer. The (layer, row) pairs
// are shared out among the OpenMP threads; each cell has exactly one writer.
inline void update_region_3d(const Field3D& c, Field3D& next, const FaceVelocities& vel, const StencilCoeffs3D& k,
                             int i0, int i1, int j0, int j1) {
    if (i0 > i1 || j0 > j1) return;
    const int nz = c.nz();
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int layer = 1; layer <= nz; ++layer) {
        for (int i = i0; i <= i1; ++i) {
            update_row_3d(c, next, vel, k, layer, i, j0, j1);
        }
    }
}
