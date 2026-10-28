// Unit test for the serial solver (src/serial/solver.cpp) on small grids:
// the correctness checks of docs/experiment_plan.md that need no parallel run.
// Built and run by CTest: ctest --test-dir build

#include "serial/solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
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

// A small basin with a developed flow: 48 x 96 cells at Re = 50 (U = 1,
// nu = 0.02, F0 = U nu (2 pi / ly)^3), after a short spin-up.
Config small_basin(int steps) {
    Config c;
    c.nx = 48;
    c.ny = 96;
    c.steps = steps;
    c.spinup = 1.0;
    c.viscosity = 0.02;
    c.wind_amplitude = 0.02 * M_PI * M_PI * M_PI;
    c.max_speed = 2.0;
    c.blob_width = 0.15;  // several cells wide even on a coarse grid
    c.blob_y = 1.0;
    c.blob_x = 0.3;
    c.check_interval = 10;
    return c;
}

}  // namespace

int main() {
    std::cout << "Pigment in the computed flow\n";
    {
        const Config c = small_basin(300);
        const SerialResult r = run_serial(c);
        check(r.summary.u_max_measured > 0.1, "the wind has spun up a flow (largest speed " +
                                                  sci(std::max(r.summary.u_max_measured, r.summary.v_max_measured)) + ")");
        check(r.summary.mass_error < 1e-12, "relative mass change < 1e-12 (got " + sci(r.summary.mass_error) + ")");
        const double lowest = *std::min_element(r.field.begin(), r.field.end());
        const double highest = *std::max_element(r.field.begin(), r.field.end());
        check(lowest >= -1e-14 && highest <= 1.0 + 1e-14, "no new extrema: the pigment stays within [0, 1]");
    }

    std::cout << "Constant field under the flow\n";
    {
        Config c = small_basin(100);
        c.initial_case = InitialCase::constant;
        c.background = 42.0;
        const SerialResult r = run_serial(c);
        double deviation = 0.0;
        for (double v : r.field) deviation = std::max(deviation, std::abs(v - 42.0));
        check(deviation <= 1e-12 * 42.0, "stays 42 to rounding (max deviation " + sci(deviation) + ")");
    }

    std::cout << "Rest stays at rest\n";
    {
        Config c = small_basin(100);
        c.wind_amplitude = 0.0;
        const SerialResult r = run_serial(c);
        check(r.summary.energy == 0.0 && r.summary.u_max_measured == 0.0 && r.summary.v_max_measured == 0.0,
              "without wind, the flow stays exactly at rest");
    }

    std::cout << "The two wind forms\n";
    {
        Config steady = small_basin(100);
        Config periodic = steady;
        periodic.wind = Wind::periodic;
        periodic.wind_alpha = 0.0;
        periodic.wind_period = 0.3;
        const SerialResult a = run_serial(steady);
        const SerialResult b = run_serial(periodic);
        check(a.field == b.field && a.summary.energy == b.summary.energy,
              "periodic wind with alpha = 0 gives exactly the steady result");
        periodic.wind_alpha = 0.5;
        const SerialResult c = run_serial(periodic);
        check(c.summary.energy != a.summary.energy, "with alpha = 0.5 the result differs");
    }

    std::cout << "Time-step check\n";
    {
        // A strong wind with a far too small speed estimate: viscosity sets dt at
        // first, but the flow outgrows it and the advection limit is exceeded.
        Config c = small_basin(100);
        c.wind_amplitude = 50.0;
        c.max_speed = 0.02;
        bool stopped = false;
        try {
            run_serial(c);
        } catch (const std::runtime_error& e) {
            stopped = std::string(e.what()).find("limit") != std::string::npos;
        }
        check(stopped, "a run with too large a dt stops with a message naming the limit");
    }

    std::cout << "Cosine decay without wind\n";
    {
        Config c;
        c.nx = 40;
        c.ny = 32;
        c.steps = 40;
        c.initial_case = InitialCase::cosine_decay;
        c.wind_amplitude = 0.0;
        c.max_speed = 0.0;
        const SerialResult r = run_serial(c);
        check(r.summary.cosine_discrete_error < 1e-12, "matches (1 - lambda dt)^n exactly");
        check(r.summary.cosine_continuous_error < 1e-3, "close to the continuous solution");
        check(r.summary.mass_error < 1e-12, "zero total mass stays zero");
    }

    // dt from cfl shrinks like dx^2 here (diffusion), so the time error O(dt)
    // and the spatial error O(dx^2) both drop about 4x when the grid is refined 2x.
    std::cout << "Cosine decay converges under grid refinement\n";
    {
        double previous = 0.0;
        bool decreasing = true;
        for (int n : {16, 32, 64}) {
            Config c;
            c.nx = n;
            c.ny = 2 * n;
            c.steps = 1;
            c.initial_case = InitialCase::cosine_decay;
            c.wind_amplitude = 0.0;
            c.max_speed = 0.0;
            c.kappa = 0.005;
            c.viscosity = 1e-9;  // irrelevant without flow; keeps its dt limit out of the way
            c.cfl = 0.5;
            const double dt = time_step(c);
            c.steps = static_cast<int>(std::lround(0.5 / dt));
            const SerialResult r = run_serial(c);
            std::cout << "    " << n << " x " << 2 * n << ": error vs continuous = " << sci(r.summary.cosine_continuous_error) << "\n";
            if (previous > 0) decreasing = decreasing && r.summary.cosine_continuous_error < previous / 3;
            previous = r.summary.cosine_continuous_error;
        }
        check(decreasing, "error drops more than 3x per refinement (second order)");
    }

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
