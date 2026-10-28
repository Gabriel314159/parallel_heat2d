// Stability limits of the time step. See docs/design.md, "Time step".
//
// With face speeds u_max, v_max, define the rates
//   a = u_max/dx + v_max/dy              (advection)
//   d = 2 kappa (1/dx^2 + 1/dy^2)        (pigment diffusion)
//   n = nu (1/dx^2 + 1/dy^2)             (viscosity)
// The time step must satisfy
//   pigment, forward Euler:      dt (a + d) <= 1
//   vorticity advection, AB3:    dt a <= 0.5
//   viscosity, AB3:              dt n <= 0.1
// AB3 is stable for centred advection up to |lambda dt| = 0.72 on the imaginary
// axis and for diffusion up to 6/11 on the negative real axis, where the
// viscous eigenvalue is -4n; the limits above keep a safety margin.

#pragma once

#include "common/config.hpp"
#include "common/grid.hpp"

// Largest dt that satisfies every limit above for face speeds u_max, v_max.
double stable_dt(const Config& c, double u_max, double v_max);

// Throws std::runtime_error naming the violated limit if dt is too large for
// the measured speeds u_max, v_max.
void check_time_step(const Config& c, double dt, double u_max, double v_max);

// Largest |U| over the east faces and largest |V| over the north faces of a block.
void max_face_speeds(const FaceVelocities& vel, double& u_max, double& v_max);
