// The vorticity equation on the cell corners, and the face velocities it gives.
// See docs/design.md, "The wind-driven flow", "Vorticity" and "One time step".
//
// Like the pigment stencil, these functions work on one row of a block at a
// time and contain no communication or threading: the serial and parallel
// solvers decide which rows to visit.

#pragma once

#include <vector>

#include "common/config.hpp"
#include "common/corners.hpp"
#include "common/grid.hpp"

// Constants of the right-hand side, computed once per run.
struct FlowCoeffs {
    double inv_12dxdy;  // 1 / (12 dx dy): the Arakawa average of three terms over 4 dx dy each
    double nu_dx2;      // nu / dx^2
    double nu_dy2;      // nu / dy^2
};

FlowCoeffs make_flow_coeffs(double nu, double dx, double dy);

// sin(2 pi y / ly) at the local corner rows i = 0..ny of a block.
std::vector<double> wind_profile(const Config& c, const Block& b);

// The time factor of the wind: F0 for the steady wind, F0 (1 + alpha sin(omega t))
// for the periodic one, with omega = 2 pi / wind_period. F = wind_strength * profile.
double wind_strength(const Config& c, double t);

// R = J(psi, zeta) + F + nu lap(zeta) at the corners (i, j0..j1).
// forcing is F on this row; extra, if not null, is an additional source per
// column (extra[j]), used by the manufactured-solution test. With
// advection = false the Jacobian is left out (the linear problem).
// psi and zeta need their halos, including the corner ghosts.
void vorticity_rhs_row(const Field& psi, const Field& zeta, Field& rhs, const FlowCoeffs& k, double forcing,
                       const double* extra, int i, int j0, int j1, bool advection = true);

// Adams–Bashforth weights for flow step number s = 0, 1, 2, ...: forward Euler
// for the first step, AB2 for the second, AB3 from then on.
struct AbWeights {
    double w0;  // for R^n
    double w1;  // for R^{n-1}
    double w2;  // for R^{n-2}
};
AbWeights ab_weights(long s);

// zeta += dt (w0 r0 + w1 r1 + w2 r2) at the corners (i, j0..j1).
void ab_update_row(Field& zeta, const Field& r0, const Field& r1, const Field& r2, const AbWeights& w, double dt,
                   int i, int j0, int j1);

// Face velocities from psi at the corners (docs/design.md, "Velocity on faces"):
//   U_e(i, j) =  (psi(i-1, j) - psi(i, j)) / dy        i = 1..ny,   j = 0..nx
//   V_n(i, j) = -(psi(i-1, j) - psi(i-1, j-1)) / dx    i = 1..ny+1, j = 1..nx
// Row i fills the east faces of row i (if i <= ny) and the north faces of row i.
// psi's ghost row 0 and ghost column 0, including the corner (0, 0), must be filled.
void faces_from_psi_row(const Field& psi, FaceVelocities& vel, double dx, double dy, int i);
void faces_from_psi(const Field& psi, FaceVelocities& vel, double dx, double dy);

// Local sums over a block's unknown corners: energy = 1/2 dx dy sum(psi zeta),
// enstrophy = 1/2 dx dy sum(zeta^2). The parallel solver adds them up across ranks.
void local_energy_enstrophy(const Field& psi, const Field& zeta, const Block& b, double dx, double dy,
                            double& energy, double& enstrophy);
