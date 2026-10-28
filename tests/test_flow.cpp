// Unit test for src/common/flow.cpp: the vorticity equation on the corners.
// See docs/experiment_plan.md, "Computed flow".
// Built and run by CTest: ctest --test-dir build

#include "common/flow.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#include "common/corners.hpp"
#include "common/grid.hpp"
#include "common/poisson.hpp"

namespace {

int failures = 0;
const double pi = M_PI;

std::string sci(double value) {
    std::ostringstream out;
    out.precision(2);
    out << std::scientific << value;
    return out.str();
}

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

Block whole(int nx, int ny) {
    Block b;
    b.nx = nx;
    b.ny = ny;
    return b;
}

// The flow on a whole 1 x 2 basin, stepped like the serial solver but with the
// forcing given per row and an optional extra source field.
struct TestFlow {
    int nx, ny;
    double lx = 1.0, ly = 2.0, dx, dy;
    Block block;
    Field psi, zeta, r0, r1, r2, source;
    FlowCoeffs k;
    CornerRange unknown;
    PoissonCG poisson;
    SerialComm comm;
    bool advection = true;
    bool use_source = false;

    TestFlow(int nx_, int ny_, double nu)
        : nx(nx_), ny(ny_), dx(1.0 / nx_), dy(2.0 / ny_), block(),
          psi(ny_, nx_), zeta(ny_, nx_), r0(ny_, nx_), r1(ny_, nx_), r2(ny_, nx_), source(ny_, nx_),
          k(make_flow_coeffs(nu, 1.0 / nx_, 2.0 / ny_)), unknown(),
          poisson(whole(nx_, ny_), 1.0 / nx_, 2.0 / ny_, 1e-13, 100000) {
        block = whole(nx_, ny_);
        unknown = unknown_corners(block);
    }

    double x(int j) const { return corner_x(j, dx); }
    double y(int i) const { return corner_y(i, dy, ly); }

    // Set zeta at the unknown corners from f(x, y), then psi from the Poisson equation.
    template <typename F>
    void set_zeta(F f) {
        for (int i = unknown.i_first; i <= unknown.i_last; ++i) {
            for (int j = unknown.j_first; j <= unknown.j_last; ++j) zeta(i, j) = f(x(j), y(i));
        }
        poisson.solve(zeta, psi, comm);
    }

    void step(long s, double dt, const std::function<double(int)>& row_forcing) {
        for (int i = unknown.i_first; i <= unknown.i_last; ++i) {
            vorticity_rhs_row(psi, zeta, r0, k, row_forcing(i), use_source ? source.row(i) : nullptr, i,
                              unknown.j_first, unknown.j_last, advection);
        }
        const AbWeights w = ab_weights(s);
        for (int i = unknown.i_first; i <= unknown.i_last; ++i) {
            ab_update_row(zeta, r0, r1, r2, w, dt, i, unknown.j_first, unknown.j_last);
        }
        poisson.solve(zeta, psi, comm);
        std::swap(r2, r1);
        std::swap(r1, r0);
    }

    void energy_enstrophy(double& e, double& z) const { local_energy_enstrophy(psi, zeta, block, dx, dy, e, z); }
};

}  // namespace

int main() {
    std::cout << "Conservation by the Arakawa Jacobian (no wind, no viscosity)\n";
    {
        // Relative drift of energy and enstrophy over the same time with dt and dt/2.
        double drift_e[2], drift_z[2];
        for (int run = 0; run < 2; ++run) {
            TestFlow f(48, 96, 0.0);
            f.set_zeta([](double x, double y) {
                return 40 * (std::sin(pi * x) * std::sin(pi * y / 2) + 0.7 * std::sin(2 * pi * x) * std::sin(3 * pi * y / 2) +
                             0.4 * std::sin(3 * pi * x) * std::sin(pi * y));
            });
            double e0, z0;
            f.energy_enstrophy(e0, z0);
            const double dt = run == 0 ? 4e-4 : 2e-4;
            const long steps = std::lround(0.4 / dt);
            for (long s = 0; s < steps; ++s) f.step(s, dt, [](int) { return 0.0; });
            double e1, z1;
            f.energy_enstrophy(e1, z1);
            drift_e[run] = std::abs(e1 - e0) / e0;
            drift_z[run] = std::abs(z1 - z0) / z0;
            std::cout << "    dt = " << dt << ": energy drift " << sci(drift_e[run]) << ", enstrophy drift "
                      << sci(drift_z[run]) << "\n";
        }
        check(drift_e[0] < 1e-4 && drift_z[0] < 1e-4, "energy and enstrophy are conserved up to the time-stepping error");
        // Halving dt must shrink the drift at least about 4x: it is time-stepping
        // error, with the forward-Euler start-up step contributing at order dt^2.
        check(drift_e[1] < drift_e[0] / 3 && drift_z[1] < drift_z[0] / 3,
              "the drift is time-stepping error: it shrinks at least like dt^2");
    }

    std::cout << "Rest stays at rest\n";
    {
        TestFlow f(16, 32, 0.01);
        for (long s = 0; s < 20; ++s) f.step(s, 1e-3, [](int) { return 0.0; });
        bool zero = true;
        for (int i = 0; i <= f.ny + 1; ++i) {
            for (int j = 0; j <= f.nx + 1; ++j) zero = zero && f.psi(i, j) == 0.0 && f.zeta(i, j) == 0.0;
        }
        check(zero, "without wind, psi and zeta stay exactly 0");
    }

    std::cout << "Linear steady state (no advection) against its series solution\n";
    {
        // Steady, without advection: nu lap(zeta) = -F and zeta = -lap(psi), with
        // psi = zeta = 0 on the walls. Expanding 1 = sum_{m odd} 4/(m pi) sin(m pi x/lx),
        //   psi = sum_{m odd} 4 F0 / (m pi nu lambda_m^2) sin(m pi x/lx) sin(2 pi y/ly),
        //   lambda_m = (m pi/lx)^2 + (2 pi/ly)^2.
        const double f0 = 1.0, nu = 0.05, lx = 1.0, ly = 2.0;
        auto series = [&](double x, double y) {
            double sum = 0.0;
            for (int m = 1; m < 400; m += 2) {
                const double lambda = (m * pi / lx) * (m * pi / lx) + (2 * pi / ly) * (2 * pi / ly);
                sum += 4 * f0 / (m * pi * nu * lambda * lambda) * std::sin(m * pi * x / lx);
            }
            return sum * std::sin(2 * pi * y / ly);
        };

        // The discrete steady state: two Poisson solves, -lap(zeta) = F / nu and -lap(psi) = zeta.
        double previous = 0.0;
        bool fourfold = true;
        Field steady32(64, 32);
        for (int nx : {32, 64}) {
            TestFlow f(nx, 2 * nx, nu);
            Field rhs(2 * nx, nx), zeta(2 * nx, nx), psi(2 * nx, nx);
            for (int i = 1; i < 2 * nx; ++i) {
                for (int j = 1; j < nx; ++j) rhs(i, j) = f0 * std::sin(2 * pi * f.y(i) / ly) / nu;
            }
            f.poisson.solve(rhs, zeta, f.comm);
            f.poisson.solve(zeta, psi, f.comm);
            double error = 0.0, size = 0.0;
            for (int i = 1; i < 2 * nx; ++i) {
                for (int j = 1; j < nx; ++j) {
                    error = std::max(error, std::abs(psi(i, j) - series(f.x(j), f.y(i))));
                    size = std::max(size, std::abs(series(f.x(j), f.y(i))));
                }
            }
            std::cout << "    " << nx << " x " << 2 * nx << ": discrete steady state vs series, relative error "
                      << sci(error / size) << "\n";
            if (previous > 0) fourfold = fourfold && previous / (error / size) > 3.5;
            previous = error / size;
            if (nx == 32) steady32 = psi;
        }
        check(fourfold, "the discrete steady state converges to the series at second order");

        // Time-stepping from rest must reach that discrete steady state.
        TestFlow f(32, 64, nu);
        f.advection = false;
        const double dt = 0.8 * 0.1 / (nu * (1 / (f.dx * f.dx) + 1 / (f.dy * f.dy)));
        const long steps = std::lround(25.0 / dt);
        for (long s = 0; s < steps; ++s) {
            f.step(s, dt, [&](int i) { return f0 * std::sin(2 * pi * f.y(i) / ly); });
        }
        double error = 0.0, size = 0.0;
        for (int i = 1; i < 64; ++i) {
            for (int j = 1; j < 32; ++j) {
                error = std::max(error, std::abs(f.psi(i, j) - steady32(i, j)));
                size = std::max(size, std::abs(steady32(i, j)));
            }
        }
        check(error / size < 1e-6, "time-stepping from rest reaches it (relative difference " + sci(error / size) + ")");

        // F0 > 0: psi > 0 in the southern half, i.e. the southern gyre turns counterclockwise.
        check(f.psi(48, 16) > 0 && f.psi(16, 16) < 0, "southern gyre counterclockwise, northern gyre clockwise");
    }

    std::cout << "Manufactured solution (every term, second order in space)\n";
    {
        // psi = a(t) S1 + b(t) S2 with two sine modes that vanish, with their
        // Laplacians, on the walls; the source makes it an exact solution.
        const double nu = 0.01, lx = 1.0, ly = 2.0;
        const double k1 = pi / lx, l1 = pi / ly;
        const double lam1 = k1 * k1 + l1 * l1, lam2 = 4 * lam1;
        auto a = [](double t) { return 1 + 0.5 * std::sin(2 * t); };
        auto da = [](double t) { return std::cos(2 * t); };
        auto b = [](double t) { return 0.5 * std::cos(3 * t); };
        auto db = [](double t) { return -1.5 * std::sin(3 * t); };
        auto zeta_exact = [&](double x, double y, double t) {
            return a(t) * lam1 * std::sin(k1 * x) * std::sin(l1 * y) +
                   b(t) * lam2 * std::sin(2 * k1 * x) * std::sin(2 * l1 * y);
        };
        auto source = [&](double x, double y, double t) {
            const double s1 = std::sin(k1 * x) * std::sin(l1 * y), s2 = std::sin(2 * k1 * x) * std::sin(2 * l1 * y);
            const double s1x = k1 * std::cos(k1 * x) * std::sin(l1 * y), s1y = l1 * std::sin(k1 * x) * std::cos(l1 * y);
            const double s2x = 2 * k1 * std::cos(2 * k1 * x) * std::sin(2 * l1 * y);
            const double s2y = 2 * l1 * std::sin(2 * k1 * x) * std::cos(2 * l1 * y);
            const double dzeta = da(t) * lam1 * s1 + db(t) * lam2 * s2;
            const double jac = a(t) * b(t) * (lam2 - lam1) * (s1x * s2y - s1y * s2x);
            const double viscous = -nu * (a(t) * lam1 * lam1 * s1 + b(t) * lam2 * lam2 * s2);
            return dzeta - jac - viscous;
        };
        const double dt = 1e-4, t_end = 0.2;
        double previous = 0.0;
        bool fourfold = true;
        for (int nx : {32, 64}) {
            TestFlow f(nx, 2 * nx, nu);
            f.use_source = true;
            f.set_zeta([&](double x, double y) { return zeta_exact(x, y, 0.0); });
            const long steps = std::lround(t_end / dt);
            for (long s = 0; s < steps; ++s) {
                for (int i = f.unknown.i_first; i <= f.unknown.i_last; ++i) {
                    for (int j = f.unknown.j_first; j <= f.unknown.j_last; ++j) f.source(i, j) = source(f.x(j), f.y(i), s * dt);
                }
                f.step(s, dt, [](int) { return 0.0; });
            }
            double error = 0.0, size = 0.0;
            for (int i = f.unknown.i_first; i <= f.unknown.i_last; ++i) {
                for (int j = f.unknown.j_first; j <= f.unknown.j_last; ++j) {
                    const double exact = zeta_exact(f.x(j), f.y(i), steps * dt);
                    error = std::max(error, std::abs(f.zeta(i, j) - exact));
                    size = std::max(size, std::abs(exact));
                }
            }
            std::cout << "    " << nx << " x " << 2 * nx << ": relative error in zeta " << sci(error / size) << "\n";
            if (previous > 0) fourfold = fourfold && previous / (error / size) > 3.5;
            previous = error / size;
        }
        check(fourfold, "error falls about 4x per refinement");
    }

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
