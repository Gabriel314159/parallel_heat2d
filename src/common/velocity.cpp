// Function bodies for common/velocity.hpp.

#include "common/velocity.hpp"

#include <cmath>

// "constexpr": a constant fixed at compile time. Being outside any function,
// it is visible to every function below in this file.
constexpr double pi = M_PI;

double streamfunction(const FlowParams& p, double x, double y, double t) {
    const double X = 2 * x / p.lx;
    const double Y = y / p.ly;
    const double f = X + (p.epsilon * std::sin(p.omega * t) * (X*X - 2*X));
    const double psi = p.amplitude * std::sin(pi * f) * std::sin(pi * Y);
    return psi;
}

double max_speed_u(const FlowParams& p) {
    const double u_max = (std::abs(p.amplitude) *  pi)/ p.ly;
    return u_max;
}

double max_speed_v(const FlowParams& p) {
    const double v_max = 2 * std::abs(p.amplitude) * pi * (1 + 2*p.epsilon) / p.lx;
    return v_max;
}

// Worked example. The east face of cell (i, j) is the vertical segment
//   x = j*dx,  from y_bottom = ly - i*dy  up to  y_top = ly - (i-1)*dy.
// U_e = (psi(top) - psi(bottom)) / dy.
double face_velocity_east(const FlowParams& p, int i, int j, double dx, double dy, double t) {
    const double x = j * dx;
    const double y_top = p.ly - (i - 1) * dy;
    const double y_bottom = p.ly - i * dy;
    return (streamfunction(p, x, y_top, t) - streamfunction(p, x, y_bottom, t)) / dy;
}

// The north face of cell (i, j) is the horizontal segment
//   y = ly - (i-1)*dy,  from x_left = (j-1)*dx  to  x_right = j*dx.
// V_n = -(psi(right) - psi(left)) / dx.
double face_velocity_north(const FlowParams& p, int i, int j, double dx, double dy, double t) {
    const double y = p.ly - (i - 1) * dy;
    const double x_left = (j - 1) * dx;
    const double x_right = j * dx;
    return -(streamfunction(p, x_right, y, t) - streamfunction(p, x_left, y, t)) / dx;
}

// The factor expressions below are written exactly as in streamfunction(),
// so corner_psi() reproduces streamfunction() bit for bit.
FlowField::FlowField(const FlowParams& p, const Block& b, double dx, double dy)
    : p_(p), b_(b), dx_(dx), dy_(dy), sx_(b.nx + 1, 0.0), sy_(b.ny + 1, 0.0) {
    for (int ic = 0; ic <= b_.ny; ++ic) {
        const double y = p_.ly - (b_.row_offset + ic) * dy_;
        const double Y = y / p_.ly;
        sy_[ic] = std::sin(pi * Y);
    }
    set_time(0.0);
}

bool FlowField::steady() const {
    return p_.amplitude == 0.0 || p_.epsilon == 0.0;
}

void FlowField::set_time(double t) {
    const double s = std::sin(p_.omega * t);
    for (int jc = 0; jc <= b_.nx; ++jc) {
        const double x = (b_.col_offset + jc) * dx_;
        const double X = 2 * x / p_.lx;
        const double f = X + (p_.epsilon * s * (X*X - 2*X));
        sx_[jc] = std::sin(pi * f);
    }
}

// East face of (i, j): corners (i-1, j) at the top and (i, j) at the bottom.
void FlowField::fill_east_row(FaceVelocities& vel, int i) const {
    double* u = vel.east_row(i);
    for (int j = 0; j <= b_.nx; ++j) {
        u[j] = (corner_psi(i - 1, j) - corner_psi(i, j)) / dy_;
    }
    if (b_.wall_left) u[0] = 0.0;
    if (b_.wall_right) u[b_.nx] = 0.0;
}

// North face of (i, j): corners (i-1, j-1) on the left and (i-1, j) on the right.
void FlowField::fill_north_row(FaceVelocities& vel, int i) const {
    double* v = vel.north_row(i);
    const bool on_wall = (i == 1 && b_.wall_top) || (i == b_.ny + 1 && b_.wall_bottom);
    for (int j = 1; j <= b_.nx; ++j) {
        v[j - 1] = on_wall ? 0.0 : -(corner_psi(i - 1, j) - corner_psi(i - 1, j - 1)) / dx_;
    }
}

void FlowField::fill_all(FaceVelocities& vel) const {
    for (int i = 1; i <= b_.ny + 1; ++i) {
        if (i <= b_.ny) fill_east_row(vel, i);
        fill_north_row(vel, i);
    }
}
