// Unit test for src/common/config.cpp. Built and run by CTest: ctest --test-dir build

#include "common/config.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

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
        check(c.nx == 256 && c.ny == 256 && c.steps == 10000, "grid and steps from defaults.json");
        check(!c.dt, "dt is automatic by default");
        check(close(time_step(c), 4.900702313605252e-4, 1e-12), "automatic dt = 4.9007e-4");
        const FlowParams p = flow_params(c);
        check(close(p.omega, 2.0 * M_PI / 10.0, 1e-15), "omega = 2 pi / period");
        check(p.amplitude == 0.1 && p.epsilon == 0.25, "flow amplitude and epsilon copied");
    }

    std::cout << "Parsing\n";
    {
        const Config c = parse({"--nx", "37", "--ny", "29", "--kappa", "0.01", "--case", "cosine_decay",
                                "--flow-amplitude", "0", "--blob-x", "0.3", "--output", "out.csv"});
        validate(c);
        check(c.nx == 37 && c.ny == 29, "--nx, --ny");
        check(c.kappa == 0.01, "--kappa");
        check(c.initial_case == InitialCase::cosine_decay, "--case cosine_decay");
        check(c.flow_amplitude == 0.0 && c.blob_x == 0.3, "--flow-amplitude, --blob-x");
        check(c.output == "out.csv", "--output");
        // Pure diffusion: dt = cfl / d with d = 2 kappa (37^2 + 29^2).
        check(close(time_step(c), 0.8 / (2 * 0.01 * (37.0 * 37.0 + 29.0 * 29.0)), 1e-12),
              "zero-flow dt = cfl / d");
    }
    {
        const Config c = parse({"--dt", "1e-4"});
        validate(c);
        check(time_step(c) == 1e-4, "explicit stable --dt is used as given");
    }

    std::cout << "Rejected runs\n";
    expect_rejected({"--dt", "1e-3"}, "unstable explicit dt");
    expect_rejected({"--dt", "-1e-4"}, "negative dt");
    expect_rejected({"--nx", "0"}, "nx = 0");
    expect_rejected({"--ly", "0"}, "ly = 0");
    expect_rejected({"--steps", "0"}, "steps = 0");
    expect_rejected({"--kappa", "-1"}, "negative kappa");
    expect_rejected({"--cfl", "0"}, "cfl = 0");
    expect_rejected({"--cfl", "1.5"}, "cfl > 1");
    expect_rejected({"--flow-epsilon", "-0.1"}, "negative epsilon");
    expect_rejected({"--flow-period", "0"}, "period = 0");
    expect_rejected({"--blob-width", "0"}, "blob width = 0");
    expect_rejected({"--case", "cosine_decay"}, "cosine_decay with the flow on");
    expect_rejected({"--case", "heated_top"}, "removed case name");
    expect_rejected({"--bogus", "3"}, "unknown flag");
    expect_rejected({"--nx"}, "flag without a value");
    expect_rejected({"--kappa", "abc"}, "value that is not a number");

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
