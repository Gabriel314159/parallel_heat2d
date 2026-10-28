// Unit test for src/common/poisson.hpp: the conjugate gradient solver of the
// Poisson equation for psi. Built and run by CTest: ctest --test-dir build

#include "common/poisson.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "common/corners.hpp"
#include "common/grid.hpp"

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

// Deterministic pseudo-random numbers in [-1, 1).
double next_random(std::uint64_t& state) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(state >> 11) / 4503599627370496.0 - 1.0;
}

Block whole(int nx, int ny) {
    Block b;
    b.nx = nx;
    b.ny = ny;
    return b;
}

// The five-point operator -(lap psi) at an unknown corner.
double minus_laplacian(const Field& psi, int i, int j, double dx, double dy) {
    return -((psi(i, j + 1) - 2 * psi(i, j) + psi(i, j - 1)) / (dx * dx) +
             (psi(i - 1, j) - 2 * psi(i, j) + psi(i + 1, j)) / (dy * dy));
}

// Solve on an nx x ny grid of a 1 x 2 basin with the right-hand side and exact
// solution given as functions of (x, y); return the largest error.
template <typename Rhs, typename Exact>
double solve_error(int nx, int ny, double tol, Rhs rhs, Exact exact) {
    const double lx = 1.0, ly = 2.0, dx = lx / nx, dy = ly / ny;
    Field zeta(ny, nx), psi(ny, nx);
    for (int i = 1; i < ny; ++i) {
        for (int j = 1; j < nx; ++j) zeta(i, j) = rhs(corner_x(j, dx), corner_y(i, dy, ly));
    }
    PoissonCG solver(whole(nx, ny), dx, dy, tol, 100000);
    solver.solve(zeta, psi, SerialComm());
    double error = 0.0;
    for (int i = 1; i < ny; ++i) {
        for (int j = 1; j < nx; ++j) {
            error = std::max(error, std::abs(psi(i, j) - exact(corner_x(j, dx), corner_y(i, dy, ly))));
        }
    }
    return error;
}

}  // namespace

int main() {
    std::cout << "Discrete eigenmode\n";
    {
        const int nx = 32, ny = 64;
        const double lx = 1.0, ly = 2.0, dx = lx / nx, dy = ly / ny;
        const double sx = std::sin(pi * dx / (2 * lx)), sy = std::sin(pi * dy / (2 * ly));
        const double lambda = 4 * sx * sx / (dx * dx) + 4 * sy * sy / (dy * dy);
        auto mode = [&](double x, double y) { return std::sin(pi * x / lx) * std::sin(pi * y / ly); };
        const double error = solve_error(nx, ny, 1e-12, [&](double x, double y) { return lambda * mode(x, y); }, mode);
        check(error < 1e-10, "sin(pi x/lx) sin(pi y/ly) is reproduced to the tolerance (error " + sci(error) + ")");
    }

    std::cout << "Stopping rule on a random right-hand side\n";
    {
        const int nx = 32, ny = 64;
        const double dx = 1.0 / nx, dy = 2.0 / ny, tol = 1e-10;
        Field zeta(ny, nx), psi(ny, nx);
        std::uint64_t state = 42;
        double bb = 0.0;
        for (int i = 1; i < ny; ++i) {
            for (int j = 1; j < nx; ++j) {
                zeta(i, j) = next_random(state);
                bb += zeta(i, j) * zeta(i, j);
            }
        }
        PoissonCG solver(whole(nx, ny), dx, dy, tol, 100000);
        const CgStats stats = solver.solve(zeta, psi, SerialComm());
        double rr = 0.0;
        for (int i = 1; i < ny; ++i) {
            for (int j = 1; j < nx; ++j) {
                const double r = zeta(i, j) - minus_laplacian(psi, i, j, dx, dy);
                rr += r * r;
            }
        }
        const double relative = std::sqrt(rr / bb);
        check(relative <= 1.5 * tol, "true relative residual " + sci(relative) + " after " +
                                         std::to_string(stats.iterations) + " iterations (tolerance " + sci(tol) + ")");
        bool walls_untouched = true;
        for (int j = 0; j <= nx; ++j) walls_untouched = walls_untouched && psi(0, j) == 0.0 && psi(ny, j) == 0.0;
        for (int i = 0; i <= ny; ++i) walls_untouched = walls_untouched && psi(i, 0) == 0.0 && psi(i, nx) == 0.0;
        check(walls_untouched, "wall corners stay 0");

        // Warm start from the solution: no iterations needed.
        const CgStats again = solver.solve(zeta, psi, SerialComm());
        check(again.iterations == 0, "starting from the solution, it stops at once");
    }

    std::cout << "Zero right-hand side\n";
    {
        Field zeta(16, 8), psi(16, 8);
        psi(3, 3) = 1.0;  // a wrong initial guess
        PoissonCG solver(whole(8, 16), 1.0 / 8, 2.0 / 16, 1e-10, 100);
        const CgStats stats = solver.solve(zeta, psi, SerialComm());
        bool zero = true;
        for (int i = 0; i <= 17; ++i) {
            for (int j = 0; j <= 9; ++j) zero = zero && psi(i, j) == 0.0;
        }
        check(zero && stats.iterations == 0, "zeta = 0 gives psi = 0 exactly");
    }

    std::cout << "Second-order convergence\n";
    {
        // psi = x sin(pi x) sin(pi y / 2) vanishes on the walls of the 1 x 2 basin.
        const double a = pi, b = pi / 2;
        auto exact = [&](double x, double y) { return x * std::sin(a * x) * std::sin(b * y); };
        auto rhs = [&](double x, double y) {
            return -2 * a * std::cos(a * x) * std::sin(b * y) + (a * a + b * b) * x * std::sin(a * x) * std::sin(b * y);
        };
        double previous = 0.0;
        bool fourfold = true;
        for (int nx : {16, 32, 64}) {
            const double error = solve_error(nx, 2 * nx, 1e-12, rhs, exact);
            std::cout << "    " << nx << " x " << 2 * nx << ": error " << sci(error) << "\n";
            if (previous > 0) fourfold = fourfold && previous / error > 3.5;
            previous = error;
        }
        check(fourfold, "error falls about 4x per refinement");
    }

    std::cout << "Too few iterations\n";
    {
        Field zeta(64, 32), psi(64, 32);
        std::uint64_t state = 5;
        for (int i = 1; i < 64; ++i) {
            for (int j = 1; j < 32; ++j) zeta(i, j) = next_random(state);
        }
        PoissonCG solver(whole(32, 64), 1.0 / 32, 2.0 / 64, 1e-12, 5);
        bool stopped = false;
        try {
            solver.solve(zeta, psi, SerialComm());
        } catch (const std::runtime_error& e) {
            stopped = std::string(e.what()).find("did not reach") != std::string::npos;
        }
        check(stopped, "a solve that cannot reach the tolerance stops with a message");
    }

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
