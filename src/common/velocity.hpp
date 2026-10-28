// The analytic double-gyre streamfunction, its speed bounds and its face velocities.
//
// The solvers no longer use this flow: they compute psi from the vorticity
// equation (common/flow.hpp, common/poisson.hpp). It is kept as an analytic
// reference for the tests of the face-velocity formulas ("Velocity on faces"
// in docs/design.md), which hold for any streamfunction.
//
// A header (.hpp) only *declares* what exists: types and function signatures.
// The function bodies live in velocity.cpp. Other files use this module by
// writing  #include "common/velocity.hpp".

// Include guard: if this header is included twice in one compilation, the
// second copy is skipped. Every header starts with this line.
#pragma once

// A struct bundles related values into one named type, like a small Python
// dataclass. Members are public and are accessed as  p.amplitude.
struct FlowParams {
    double amplitude;  // A
    double epsilon;    // oscillation strength, 0 for the steady gyre
    double omega;      // angular frequency, 2*pi / period
    double lx;         // domain width
    double ly;         // domain height
};

// "const FlowParams& p" means: pass p by reference (no copy is made) and
// promise not to modify it. This is the usual way to pass structs.

// psi(x, y, t) in physical coordinates.
double streamfunction(const FlowParams& p, double x, double y, double t);

// Analytic bounds on |u| and |v| over the whole domain and all times.
double max_speed_u(const FlowParams& p);
double max_speed_v(const FlowParams& p);

// Normal velocity on the east face of cell (i, j) and on the north face of
// cell (i, j), at time t. i and j are GLOBAL 1-based indices (rows top-down),
// so every rank computes the same value for a shared face.
double face_velocity_east(const FlowParams& p, int i, int j, double dx, double dy, double t);
double face_velocity_north(const FlowParams& p, int i, int j, double dx, double dy, double t);
