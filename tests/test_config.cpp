// Unit test for src/common/config.cpp and src/common/stability.cpp.
// Built and run by CTest: ctest --test-dir build

#include "common/config.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "common/stability.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

bool close(double a, double b, double rtol) {
    return std::abs(a - b) <= rtol * std::abs(b);
}

// Pretend these words were typed after the program name, then parse them.
Config parse(std::vector<std::string> words) {
    words.insert(words.begin(), "heat_serial");
    std::vector<char*> argv;
    for (std::string& w : words) argv.push_back(w.data());
    return parse_args(static_cast<int>(argv.size()), argv.data());
}

// A run that must be rejected: parsing, validation or the time step throws.
void expect_rejected(const std::vector<std::string>& words, const std::string& label) {
    try {
        const Config c = parse(words);
        validate(c);
        time_step(c);
    } catch (const std::exception& e) {
        check(true, label + "   [" + e.what() + "]");
        return;
    }
    check(false, label + "   (accepted, but should be rejected)");
}

}  // namespace

int main() {
    std::cout << "Defaults\n";
    {
        const Config c = parse({});
        validate(c);
        check(c.nx == 256 && c.ny == 512 && c.lx == 1.0 && c.ly == 2.0, "1 x 2 basin on 256 x 512 cells");
        check(c.viscosity == 1e-3 && c.wind_amplitude == 0.1 && c.spinup == 100.0 && c.max_speed == 2.0,
              "flow defaults nu = 1e-3, F0 = 0.1, spin-up 100, speed estimate 2");
        check(!c.dt, "dt is automatic by default");

        // The advection limit dt a <= 0.5 decides for these values.
        const double dx = 1.0 / 256, dy = 2.0 / 512;
        const double a = c.max_speed / dx + c.max_speed / dy;
        check(close(time_step(c), c.cfl * 0.5 / a, 1e-12), "automatic dt = cfl x 0.5 / a");
    }

    std::cout << "Stability limits\n";
    {
        Config c;
        c.nx = 64;
        c.ny = 128;
        const double dx = 1.0 / 64, dy = 2.0 / 128, inverse_squares = 1 / (dx * dx) + 1 / (dy * dy);
        // No flow: pigment diffusion or viscosity decides.
        const double expected = std::min(1.0 / (2 * c.kappa * inverse_squares), 0.1 / (c.viscosity * inverse_squares));
        check(close(stable_dt(c, 0.0, 0.0), expected, 1e-14), "without flow, the smallest of the other limits");
        bool rejected = false;
        try {
            check_time_step(c, 2 * stable_dt(c, 30.0, 30.0), 30.0, 30.0);
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check(rejected, "check_time_step rejects twice the stable dt");
        bool accepted = true;
        try {
            check_time_step(c, stable_dt(c, 30.0, 30.0), 30.0, 30.0);
        } catch (const std::runtime_error&) {
            accepted = false;
        }
        check(accepted, "check_time_step accepts the stable dt itself");
    }

    std::cout << "Parsing\n";
    {
        const Config c = parse({"--nx", "37", "--ny", "29", "--kappa", "0.01", "--case", "cosine_decay",
                                "--wind-amplitude", "0", "--blob-x", "0.3", "--output", "out.csv",
                                "--viscosity", "0.02", "--wind", "periodic", "--wind-alpha", "0.5",
                                "--wind-period", "3", "--spinup", "1.5", "--max-speed", "10", "--check-interval", "7"});
        validate(c);
        check(c.nx == 37 && c.ny == 29, "--nx, --ny");
        check(c.kappa == 0.01, "--kappa");
        check(c.initial_case == InitialCase::cosine_decay, "--case cosine_decay");
        check(c.wind_amplitude == 0.0 && c.blob_x == 0.3, "--wind-amplitude, --blob-x");
        check(c.output == "out.csv", "--output");
        check(c.viscosity == 0.02, "--viscosity");
        check(c.wind == Wind::periodic && c.wind_alpha == 0.5 && c.wind_period == 3, "--wind periodic, --wind-alpha, --wind-period");
        check(c.spinup == 1.5 && c.max_speed == 10 && c.check_interval == 7, "--spinup, --max-speed, --check-interval");
        check(spinup_steps(c, 0.1) == 15, "spin-up steps = spinup / dt");
    }
    {
        // --threads is optional: empty means "OMP_NUM_THREADS, or one thread".
        check(!parse({}).threads, "threads is unset by default");
        const Config c = parse({"--threads", "4"});
        validate(c);
        check(c.threads && *c.threads == 4, "--threads 4");
    }
    {
        Config c = parse({"--max-speed", "10"});
        const double limit = stable_dt(c, 10, 10);
        c = parse({"--max-speed", "10", "--dt", std::to_string(limit / 2)});
        validate(c);
        check(c.dt && time_step(c) == *c.dt && *c.dt < limit, "explicit stable --dt is used as given");
    }

    std::cout << "Rejected runs\n";
    expect_rejected({"--dt", "1"}, "unstable explicit dt");
    expect_rejected({"--dt", "-1e-4"}, "negative dt");
    expect_rejected({"--nx", "1"}, "nx = 1 (no interior corner)");
    expect_rejected({"--ly", "0"}, "ly = 0");
    expect_rejected({"--steps", "0"}, "steps = 0");
    expect_rejected({"--kappa", "-1"}, "negative kappa");
    expect_rejected({"--cfl", "0"}, "cfl = 0");
    expect_rejected({"--cfl", "1.5"}, "cfl > 1");
    expect_rejected({"--viscosity", "0"}, "zero viscosity (the only dissipation)");
    expect_rejected({"--beta", "1"}, "removed beta option");
    expect_rejected({"--wind", "gusty"}, "unknown wind");
    expect_rejected({"--wind-period", "0"}, "wind period = 0");
    expect_rejected({"--spinup", "-1"}, "negative spin-up");
    expect_rejected({"--max-speed", "-1"}, "negative speed estimate");
    expect_rejected({"--blob-width", "0"}, "blob width = 0");
    expect_rejected({"--case", "cosine_decay"}, "cosine_decay with the wind on");
    expect_rejected({"--case", "heated_top"}, "removed case name");
    expect_rejected({"--flow-amplitude", "0.1"}, "removed double-gyre option");
    expect_rejected({"--threads", "0"}, "zero threads");
    expect_rejected({"--threads", "-2"}, "negative threads");
    expect_rejected({"--threads", "many"}, "threads that is not a number");
    expect_rejected({"--bogus", "3"}, "unknown flag");
    expect_rejected({"--nx"}, "flag without a value");
    expect_rejected({"--kappa", "abc"}, "value that is not a number");

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
