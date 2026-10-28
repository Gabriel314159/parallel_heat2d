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
