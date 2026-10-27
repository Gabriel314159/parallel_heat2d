// Function bodies for common/config.hpp.

#include "common/config.hpp"

#include <algorithm>  // std::max
#include <cmath>
#include <stdexcept>  // std::invalid_argument
#include <string>

constexpr double pi = M_PI;

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
        "  --nx N, --ny N          cells along x and y (default 256 x 256)\n"
        "  --lx L, --ly L          domain size (default 1 x 1)\n"
        "  --steps N               number of time steps (default 10000)\n"
        "  --kappa K               diffusivity (default 0.005)\n"
        "  --cfl C                 safety fraction of the stable dt, in (0, 1] (default 0.8)\n"
        "  --dt DT                 explicit time step instead of --cfl\n"
        "\n"
        "Double-gyre flow\n"
        "  --flow-amplitude A      flow strength; 0 switches the flow off (default 0.1)\n"
        "  --flow-epsilon E        oscillation strength; 0 = steady (default 0.25)\n"
        "  --flow-period T         oscillation period (default 10)\n"
        "\n"
        "Initial condition\n"
        "  --case NAME             pigment_drop | constant | cosine_decay (default pigment_drop)\n"
        "  --background C0         background value (default 0)\n"
        "  --amplitude A0          blob / cosine amplitude (default 1)\n"
        "  --blob-x X, --blob-y Y  blob centre (default 0.5, 0.5)\n"
        "  --blob-width SIGMA      blob width (default 0.05)\n"
        "\n"
        "Output\n"
        "  --output PATH           write the final field: PATH.csv as text, PATH.npy as binary;\n"
        "                          the run summary goes next to it as .json\n"
        "  --snapshot-interval N   write a snapshot every N steps (default 0 = off)\n"
        "  --snapshot-dir DIR      snapshot folder (default outputs/snapshots)\n"
        "  --snapshot-coarsen K    average snapshots over K x K cells to keep them small (default 1)\n"
        "  -h, --help              show this help\n"
        "\n"
        "Parallel execution (heat_parallel only; threads per rank come from OMP_NUM_THREADS)\n"
        "  --decomp rows|2d        split rows only, or a 2D process grid (default rows)\n"
        "  --exchange blocking|overlap\n"
        "                          blocking halo exchange, or nonblocking with computation\n"
        "                          of the interior while messages travel (default blocking)\n"
        "  --py N, --px N          process grid rows and columns for --decomp 2d\n"
        "                          (default: chosen automatically)\n";
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
        else if (flag == "--flow-amplitude") c.flow_amplitude = std::stod(value_after(k, argc, argv));
        else if (flag == "--flow-epsilon") c.flow_epsilon = std::stod(value_after(k, argc, argv));
        else if (flag == "--flow-period") c.flow_period = std::stod(value_after(k, argc, argv));
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
        else throw std::invalid_argument("unknown option: " + flag);
    }
    return c;
}

void validate(const Config& c) {
    // Worked example: one rule, one clear message.
    if (c.nx <= 0 || c.ny <= 0) {
        throw std::invalid_argument("nx and ny must be positive");
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
    if (c.flow_epsilon < 0) {
        throw std::invalid_argument("flow_epsilon must be non-negative");
    }
    if (c.flow_period <= 0) {
        throw std::invalid_argument("flow_period must be positive");
    }
    if (c.blob_width <= 0) {
        throw std::invalid_argument("blob_width must be positive");
    }
    // cosine_decay is the zero-flow analytic reference.
    if (c.initial_case == InitialCase::cosine_decay && c.flow_amplitude != 0) {
        throw std::invalid_argument("cosine_decay requires flow_amplitude = 0");
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
    if (!c.output.empty()) {
        const std::string::size_type dot = c.output.rfind('.');
        const std::string extension = dot == std::string::npos ? "" : c.output.substr(dot);
        if (extension != ".csv" && extension != ".npy") {
            throw std::invalid_argument("--output must end in .csv or .npy");
        }
    }
}

FlowParams flow_params(const Config& c) {
    const double omega = 2 * pi / c.flow_period;
    return FlowParams{c.flow_amplitude, c.flow_epsilon, omega, c.lx, c.ly};
}

double time_step(const Config& c) {
    // docs/design.md, "Time step": dt (a + d) <= 1 with
    //   a = u_max/dx + v_max/dy,  d = 2 kappa (1/dx^2 + 1/dy^2).
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;

    // max_speed_u takes a FlowParams, so build one first.
    const FlowParams p = flow_params(c);
    const double u_max = max_speed_u(p);
    const double v_max = max_speed_v(p);

    const double a = u_max / dx + v_max / dy;
    const double d = 2 * c.kappa * (1 / (dx * dx) + 1 / (dy * dy));

    if (c.dt) {
        if (*c.dt * (a + d) > 1) {
            throw std::invalid_argument("dt = " + std::to_string(*c.dt) +
                                        " is unstable; the limit is " + std::to_string(1 / (a + d)));
        }
        return *c.dt;
    }
    return c.cfl / (a + d);
}

double peclet(const Config& c) {
    const FlowParams p = flow_params(c);
    const double speed = std::max(max_speed_u(p), max_speed_v(p));
    const double length = std::max(c.lx, c.ly);
    return speed * length / c.kappa;
}

double peclet_grid(const Config& c) {
    const FlowParams p = flow_params(c);
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;
    return std::max(max_speed_u(p) * dx, max_speed_v(p) * dy) / c.kappa;
}
