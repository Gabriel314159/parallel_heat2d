// Double-gyre flow, including streamfunction, speed bounds, and face velocities.
// See docs/design.md, sections "Model" and "Velocity on faces".
//
// A header (.hpp) only *declares* what exists: types and function signatures.
// The function bodies live in velocity.cpp. Other files use this module by
// writing  #include "common/velocity.hpp".

// Include guard: if this header is included twice in one compilation, the
// second copy is skipped. Every header starts with this line.
#pragma once

#include <vector>

#include "common/grid.hpp"

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

// Fills all face velocities of a block, fast (see src/README.md, "Velocity storage").
//
// psi is separable:  psi = A * sx(x, t) * sy(y)  with  sx = sin(pi f(X, t)),
// sy = sin(pi Y). So instead of two sin() calls per face, sx is evaluated once
// per corner column (per time) and sy once per corner row, and each face
// velocity is a difference of two corner values. The arithmetic is the same
// as in face_velocity_east / face_velocity_north, so the results agree with
// them exactly, except that faces on a physical wall are set to exactly 0
// (sin(2 pi) is about -2.4e-16, not 0, in floating point).
//
// All positions come from global indices, so every block of a decomposition
// computes bit-identical values for the faces it shares with a neighbour.
class FlowField {
public:
    FlowField(const FlowParams& p, const Block& b, double dx, double dy);

    // True when the velocity does not change in time (A = 0 or epsilon = 0):
    // the faces then only need to be filled once.
    bool steady() const;

    // Recompute the corner-column factors sx for time t.
    void set_time(double t);

    // east(i, 0..nx) for one row i = 1..ny, at the time of the last set_time.
    void fill_east_row(FaceVelocities& vel, int i) const;

    // north(i, 1..nx) for one row i = 1..ny+1, at the time of the last set_time.
    void fill_north_row(FaceVelocities& vel, int i) const;

    // Every face of the block (set_time first).
    void fill_all(FaceVelocities& vel) const;

private:
    // psi at local corner (ic, jc): corner row ic = 0..ny, corner column jc = 0..nx.
    double corner_psi(int ic, int jc) const { return p_.amplitude * sx_[jc] * sy_[ic]; }

    FlowParams p_;
    Block b_;
    double dx_;
    double dy_;
    std::vector<double> sx_;  // sin(pi f(X, t)) per corner column
    std::vector<double> sy_;  // sin(pi Y) per corner row
};
