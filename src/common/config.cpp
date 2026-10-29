// Function bodies for common/config.hpp.

#include "common/config.hpp"

#include <algorithm>  // std::max
#include <cmath>
#include <stdexcept>  // std::invalid_argument
#include <string>

#include "common/stability.hpp"

// An unnamed namespace makes these helpers private to this file: other files
// cannot see or call them. Use it for small internal helper functions.
namespace {

// Return the text after the flag at position k and advance k past it.
// "int& k" is a reference: changing k here changes the caller's variable.
std::string value_after(int& k, int argc, char** argv) {
    const std::string flag = argv[k];
    if (k + 1 >= argc) {
        throw std::invalid_argument("missing value after " + flag);
    }
    ++k;
    return argv[k];
}

InitialCase parse_case(const std::string& name) {
    if (name == "pigment_drop") return InitialCase::pigment_drop;
    if (name == "constant") return InitialCase::constant;
    if (name == "cosine_decay") return InitialCase::cosine_decay;
    throw std::invalid_argument("unknown case: " + name);
}

Wind parse_wind(const std::string& name) {
    if (name == "steady") return Wind::steady;
    if (name == "periodic") return Wind::periodic;
    throw std::invalid_argument("unknown wind: " + name + " (use steady or periodic)");
}

Decomp parse_decomp(const std::string& name) {
    if (name == "rows") return Decomp::rows;
    if (name == "2d") return Decomp::two_d;
    throw std::invalid_argument("unknown decomposition: " + name + " (use rows or 2d)");
}

Exchange parse_exchange(const std::string& name) {
    if (name == "blocking") return Exchange::blocking;
    if (name == "overlap") return Exchange::overlap;
    throw std::invalid_argument("unknown exchange: " + name + " (use blocking or overlap)");
}

}  // namespace

std::string usage(const std::string& program) {
    return "Usage: " + program + " [options]\n"
        "\n"
        "Grid and time stepping\n"
        "  --nx N, --ny N          cells along x and y (default 256 x 512)\n"
        "  --lx L, --ly L          basin size (default 1 x 2)\n"
        "  --steps N               time steps after the pigment is released (default 50000)\n"
        "  --kappa K               pigment diffusivity (default 1e-4)\n"
        "  --cfl C                 safety fraction of the stable dt, in (0, 1] (default 0.8)\n"
        "  --dt DT                 explicit time step instead of --cfl and --max-speed\n"
        "  --max-speed U           estimated largest flow speed; sets the stable dt (default 2)\n"
        "  --check-interval N      check the actual speeds every N steps; 0 = never (default 100)\n"
        "  --poisson-tol TOL       Poisson solver (CG): stop at residual <= TOL x |rhs| (default 1e-10)\n"
        "  --poisson-max-iter N    Poisson solver: give up after N iterations (default 20000)\n"
        "\n"
        "Wind-driven flow\n"
        "  --viscosity NU          viscosity, > 0 (default 1e-3)\n"
        "  --wind steady|periodic  F0 sin(2 pi y/ly), or F0 (1 + alpha sin(2 pi t/T)) sin(2 pi y/ly)\n"
        "  --wind-amplitude F0     wind forcing amplitude; 0 = no wind (default 0.1)\n"
        "  --wind-alpha A          oscillation strength of the periodic wind (default 0)\n"
        "  --wind-period T         period of the periodic wind (default 1)\n"
        "  --spinup T              flow-only time before the pigment is released (default 100)\n"
        "\n"
        "Initial pigment\n"
        "  --case NAME             pigment_drop | constant | cosine_decay (default pigment_drop)\n"
        "  --background C0         background value (default 0)\n"
        "  --amplitude A0          blob / cosine amplitude (default 1)\n"
        "  --blob-x X, --blob-y Y  blob centre (default 0.2, 1.0)\n"
        "  --blob-width SIGMA      blob width (default 0.05)\n"
        "\n"
        "Output\n"
        "  --output PATH           write the final pigment field: PATH.csv as text, PATH.npy as binary;\n"
        "                          the run summary goes next to it as .json\n"
        "  --snapshot-interval N   write snapshots (pigment, psi, vorticity) every N steps (default 0 = off)\n"
        "  --snapshot-dir DIR      snapshot folder (default outputs/snapshots)\n"
        "  --snapshot-coarsen K    shrink snapshots by K in each direction (default 1)\n"
        "  -h, --help              show this help\n"
        "\n"
        "Parallel execution (heat_parallel only)\n"
        "  --decomp rows|2d        split rows only, or a 2D process grid (default rows)\n"
        "  --exchange blocking|overlap\n"
        "                          blocking halo exchanges (pigment and CG), or nonblocking with\n"
        "                          computation of the interior while messages travel (default blocking)\n"
        "  --py N, --px N          process grid rows and columns for --decomp 2d\n"
        "                          (default: chosen automatically)\n"
        "  --threads N             OpenMP threads per rank; overrides OMP_NUM_THREADS\n"
        "                          (default: OMP_NUM_THREADS, or 1 if it is not set)\n";
}

std::string wind_name(Wind w) {
    return w == Wind::steady ? "steady" : "periodic";
}

std::string decomp_name(Decomp d) {
    return d == Decomp::rows ? "rows" : "2d";
}

std::string exchange_name(Exchange e) {
    return e == Exchange::blocking ? "blocking" : "overlap";
}

bool has_parallel_options(const Config& c) {
    return c.decomposition != Decomp::rows || c.exchange != Exchange::blocking || c.py != 0 || c.px != 0;
}

bool wants_help(int argc, char** argv) {
    for (int k = 1; k < argc; ++k) {
        const std::string word = argv[k];
        if (word == "-h" || word == "--help") return true;
    }
    return false;
}

std::string case_name(InitialCase c) {
    switch (c) {
        case InitialCase::pigment_drop: return "pigment_drop";
        case InitialCase::constant: return "constant";
        case InitialCase::cosine_decay: return "cosine_decay";
    }
    return "unknown";
}

// argv is the list of words typed on the command line; argv[0] is the program
// name. std::stoi / std::stod convert text to int / double and throw if the
// text is not a number.
Config parse_args(int argc, char** argv) {
    Config c;
    for (int k = 1; k < argc; ++k) {
        const std::string flag = argv[k];
        if (flag == "--nx") c.nx = std::stoi(value_after(k, argc, argv));
        else if (flag == "--ny") c.ny = std::stoi(value_after(k, argc, argv));
        else if (flag == "--lx") c.lx = std::stod(value_after(k, argc, argv));
        else if (flag == "--ly") c.ly = std::stod(value_after(k, argc, argv));
        else if (flag == "--steps") c.steps = std::stoi(value_after(k, argc, argv));
        else if (flag == "--kappa") c.kappa = std::stod(value_after(k, argc, argv));
        else if (flag == "--cfl") c.cfl = std::stod(value_after(k, argc, argv));
        else if (flag == "--dt") c.dt = std::stod(value_after(k, argc, argv));
        else if (flag == "--max-speed") c.max_speed = std::stod(value_after(k, argc, argv));
        else if (flag == "--check-interval") c.check_interval = std::stoi(value_after(k, argc, argv));
        else if (flag == "--poisson-tol") c.poisson_tol = std::stod(value_after(k, argc, argv));
        else if (flag == "--poisson-max-iter") c.poisson_max_iter = std::stoi(value_after(k, argc, argv));
        else if (flag == "--viscosity") c.viscosity = std::stod(value_after(k, argc, argv));
        else if (flag == "--wind") c.wind = parse_wind(value_after(k, argc, argv));
        else if (flag == "--wind-amplitude") c.wind_amplitude = std::stod(value_after(k, argc, argv));
        else if (flag == "--wind-alpha") c.wind_alpha = std::stod(value_after(k, argc, argv));
        else if (flag == "--wind-period") c.wind_period = std::stod(value_after(k, argc, argv));
        else if (flag == "--spinup") c.spinup = std::stod(value_after(k, argc, argv));
        else if (flag == "--case") c.initial_case = parse_case(value_after(k, argc, argv));
        else if (flag == "--background") c.background = std::stod(value_after(k, argc, argv));
        else if (flag == "--amplitude") c.amplitude = std::stod(value_after(k, argc, argv));
        else if (flag == "--blob-x") c.blob_x = std::stod(value_after(k, argc, argv));
        else if (flag == "--blob-y") c.blob_y = std::stod(value_after(k, argc, argv));
        else if (flag == "--blob-width") c.blob_width = std::stod(value_after(k, argc, argv));
        else if (flag == "--output") c.output = value_after(k, argc, argv);
        else if (flag == "--snapshot-interval") c.snapshot_interval = std::stoi(value_after(k, argc, argv));
        else if (flag == "--snapshot-dir") c.snapshot_dir = value_after(k, argc, argv);
        else if (flag == "--snapshot-coarsen") c.snapshot_coarsen = std::stoi(value_after(k, argc, argv));
        else if (flag == "--decomp") c.decomposition = parse_decomp(value_after(k, argc, argv));
        else if (flag == "--exchange") c.exchange = parse_exchange(value_after(k, argc, argv));
        else if (flag == "--py") c.py = std::stoi(value_after(k, argc, argv));
        else if (flag == "--px") c.px = std::stoi(value_after(k, argc, argv));
        else if (flag == "--threads") c.threads = std::stoi(value_after(k, argc, argv));
        else throw std::invalid_argument("unknown option: " + flag);
    }
    return c;
}

void validate(const Config& c) {
    // At least one interior corner is needed for the flow.
    if (c.nx < 2 || c.ny < 2) {
        throw std::invalid_argument("nx and ny must be at least 2");
    }
    if (c.lx <= 0 || c.ly <= 0) {
        throw std::invalid_argument("lx and ly must be positive");
    }
    if (c.steps <= 0) {
        throw std::invalid_argument("steps must be positive");
    }
    if (c.kappa <= 0) {
        throw std::invalid_argument("kappa must be positive");
    }
    if (c.cfl <= 0 || c.cfl > 1) {
        throw std::invalid_argument("cfl must be in (0, 1]");
    }
    // c.dt is true only when the optional holds a value; *c.dt reads it.
    // "&&" stops early, so *c.dt is never read when dt is empty.
    if (c.dt && *c.dt <= 0) {
        throw std::invalid_argument("dt must be positive");
    }
    if (c.max_speed < 0) {
        throw std::invalid_argument("max_speed must be non-negative");
    }
    if (c.check_interval < 0) {
        throw std::invalid_argument("check_interval must be non-negative");
    }
    if (c.poisson_tol <= 0 || c.poisson_tol >= 1) {
        throw std::invalid_argument("poisson_tol must be in (0, 1)");
    }
    if (c.poisson_max_iter < 1) {
        throw std::invalid_argument("poisson_max_iter must be at least 1");
    }
    // Viscosity is the only dissipation: without it the wind's energy has nowhere to go.
    if (c.viscosity <= 0) {
        throw std::invalid_argument("viscosity must be positive (it is the only dissipation)");
    }
    if (c.wind_period <= 0) {
        throw std::invalid_argument("wind_period must be positive");
    }
    if (c.spinup < 0) {
        throw std::invalid_argument("spinup must be non-negative");
    }
    if (c.blob_width <= 0) {
        throw std::invalid_argument("blob_width must be positive");
    }
    // cosine_decay is the zero-flow analytic reference: no wind keeps the flow at rest.
    if (c.initial_case == InitialCase::cosine_decay && c.wind_amplitude != 0) {
        throw std::invalid_argument("cosine_decay requires wind_amplitude = 0");
    }
    if (c.snapshot_interval < 0) {
        throw std::invalid_argument("snapshot_interval must be non-negative");
    }
    if (c.snapshot_interval > 0 && c.snapshot_dir.empty()) {
        throw std::invalid_argument("snapshot_dir must not be empty");
    }
    if (c.snapshot_coarsen < 1) {
        throw std::invalid_argument("snapshot_coarsen must be at least 1");
    }
    if (c.snapshot_interval > 0 && (c.nx % c.snapshot_coarsen != 0 || c.ny % c.snapshot_coarsen != 0)) {
        throw std::invalid_argument("snapshot_coarsen must divide nx and ny");
    }
    if (c.py < 0 || c.px < 0 || (c.py > 0) != (c.px > 0)) {
        throw std::invalid_argument("--py and --px must be given together, as positive numbers");
    }
    if (c.py > 0 && c.decomposition != Decomp::two_d) {
        throw std::invalid_argument("--py and --px apply to --decomp 2d only");
    }
    if (c.threads && *c.threads < 1) {
        throw std::invalid_argument("--threads must be at least 1");
    }
    if (!c.output.empty()) {
        const std::string::size_type dot = c.output.rfind('.');
        const std::string extension = dot == std::string::npos ? "" : c.output.substr(dot);
        if (extension != ".csv" && extension != ".npy") {
            throw std::invalid_argument("--output must end in .csv or .npy");
        }
    }
}

double time_step(const Config& c) {
    // docs/design.md, "Time step": the step must satisfy every stability limit
    // for the largest speed we expect, max_speed.
    const double largest = stable_dt(c, c.max_speed, c.max_speed);
    if (c.dt) {
        if (*c.dt > largest) {
            throw std::invalid_argument("dt = " + std::to_string(*c.dt) + " is unstable for --max-speed " +
                                        std::to_string(c.max_speed) + "; the limit is " + std::to_string(largest));
        }
        return *c.dt;
    }
    return c.cfl * largest;
}

int spinup_steps(const Config& c, double dt) {
    return static_cast<int>(std::lround(c.spinup / dt));
}

double peclet(const Config& c) {
    return c.max_speed * std::max(c.lx, c.ly) / c.kappa;
}

double peclet_grid(const Config& c, double speed) {
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;
    return speed * std::max(dx, dy) / c.kappa;
}
