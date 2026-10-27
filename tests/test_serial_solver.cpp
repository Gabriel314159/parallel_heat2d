// Unit test for the serial solver (src/serial/solver.cpp) on small grids:
// the correctness checks of docs/experiment_plan.md that need no parallel run.
// Built and run by CTest: ctest --test-dir build

#include "serial/solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

namespace {

int failures = 0;

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

Config small_grid(int nx, int ny, int steps) {
    Config c;
    c.nx = nx;
    c.ny = ny;
    c.steps = steps;
    c.blob_width = 0.15;  // several cells wide even on a coarse grid
    return c;
}

}  // namespace

int main() {
    std::cout << "Mass conservation with the time-dependent double gyre\n";
    {
        const Config c = small_grid(37, 29, 200);
        const SerialResult r = run_serial(c);
        check(r.summary.mass_error < 1e-12, "relative mass change < 1e-12 (got " + sci(r.summary.mass_error) + ")");
        const double lowest = *std::min_element(r.field.begin(), r.field.end());
        const double highest = *std::max_element(r.field.begin(), r.field.end());
        check(lowest >= -1e-14 && highest <= 1.0 + 1e-14, "no new extrema: field stays within [0, 1]");
    }

    std::cout << "Mass conservation with the steady double gyre (velocities filled once)\n";
    {
        Config c = small_grid(37, 29, 200);
        c.flow_epsilon = 0.0;
        const SerialResult r = run_serial(c);
        check(r.summary.mass_error < 1e-12, "relative mass change < 1e-12");
    }

    std::cout << "Constant field under the flow\n";
    {
        Config c = small_grid(31, 23, 100);
        c.initial_case = InitialCase::constant;
        c.background = 42.0;
        const SerialResult r = run_serial(c);
        double deviation = 0.0;
        for (double v : r.field) deviation = std::max(deviation, std::abs(v - 42.0));
        check(deviation <= 1e-12 * 42.0, "stays 42 to rounding (max deviation " + sci(deviation) + ")");
    }

    std::cout << "Cosine decay without flow\n";
    {
        Config c = small_grid(40, 32, 40);
        c.initial_case = InitialCase::cosine_decay;
        c.flow_amplitude = 0.0;
        const SerialResult r = run_serial(c);
        check(r.summary.cosine_discrete_error < 1e-12, "matches (1 - lambda dt)^n exactly");
        check(r.summary.cosine_continuous_error < 1e-3, "close to the continuous solution");
        check(r.summary.mass_error < 1e-12, "zero total mass stays zero");
    }

    // dt from cfl shrinks like dx^2, so the time error O(dt) and the spatial
    // error O(dx^2) both drop about 4x when the grid is refined 2x.
    std::cout << "Cosine decay converges under grid refinement\n";
    {
        // Common final time; dt from cfl so it shrinks with dx^2.
        double previous = 0.0;
        bool decreasing = true;
        for (int n : {16, 32, 64}) {
            Config c = small_grid(n, n, 1);
            c.initial_case = InitialCase::cosine_decay;
            c.flow_amplitude = 0.0;
            c.cfl = 0.5;
            const double dt = time_step(c);
            c.steps = static_cast<int>(std::lround(0.5 / dt));
            const SerialResult r = run_serial(c);
            std::cout << "    n = " << n << ": error vs continuous = " << r.summary.cosine_continuous_error << "\n";
            if (previous > 0) decreasing = decreasing && r.summary.cosine_continuous_error < previous / 3;
            previous = r.summary.cosine_continuous_error;
        }
        check(decreasing, "error drops more than 3x per refinement (second order)");
    }

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
