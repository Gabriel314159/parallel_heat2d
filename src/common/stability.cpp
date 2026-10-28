// Function bodies for common/stability.hpp.

#include "common/stability.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

struct Rates {
    double advection;
    double diffusion;
    double viscosity;
};

Rates rates(const Config& c, double u_max, double v_max) {
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;
    const double inverse_squares = 1 / (dx * dx) + 1 / (dy * dy);
    return Rates{u_max / dx + v_max / dy, 2 * c.kappa * inverse_squares, c.viscosity * inverse_squares};
}

// limit / rate, or infinity when the rate is zero.
double allowed(double limit, double rate) {
    return rate > 0 ? limit / rate : std::numeric_limits<double>::infinity();
}

}  // namespace

double stable_dt(const Config& c, double u_max, double v_max) {
    const Rates r = rates(c, u_max, v_max);
    return std::min({allowed(1.0, r.advection + r.diffusion), allowed(0.5, r.advection), allowed(0.1, r.viscosity)});
}

void check_time_step(const Config& c, double dt, double u_max, double v_max) {
    const Rates r = rates(c, u_max, v_max);
    std::ostringstream problem;
    if (dt * (r.advection + r.diffusion) > 1.0) {
        problem << "pigment limit dt (a + d) = " << dt * (r.advection + r.diffusion) << " > 1";
    } else if (dt * r.advection > 0.5) {
        problem << "vorticity advection limit dt a = " << dt * r.advection << " > 0.5";
    } else if (dt * r.viscosity > 0.1) {
        problem << "viscosity limit dt nu (1/dx^2 + 1/dy^2) = " << dt * r.viscosity << " > 0.1";
    } else {
        return;
    }
    problem << " at measured speeds u_max = " << u_max << ", v_max = " << v_max
            << ". Rerun with a larger --max-speed or a smaller --dt.";
    throw std::runtime_error(problem.str());
}

void max_face_speeds(const FaceVelocities& vel, double& u_max, double& v_max) {
    u_max = 0.0;
    v_max = 0.0;
    for (int i = 1; i <= vel.ny(); ++i) {
        const double* u = vel.east_row(i);
        for (int j = 0; j <= vel.nx(); ++j) u_max = std::max(u_max, std::abs(u[j]));
    }
    for (int i = 1; i <= vel.ny() + 1; ++i) {
        const double* v = vel.north_row(i);
        for (int j = 0; j < vel.nx(); ++j) v_max = std::max(v_max, std::abs(v[j]));
    }
}
