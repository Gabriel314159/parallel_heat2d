// Run configuration: command-line parsing, validation and the global time step.
// Flags and defaults follow config/README.md and config/defaults.json.

#pragma once

#include <optional>
#include <string>

#include "common/velocity.hpp"

// "enum class" is a type with a fixed set of named values. A variable of this
// type can only hold one of them, so a typo like "pigment_dorp" cannot sneak
// through as it could with a plain string.
enum class InitialCase { pigment_drop, constant, cosine_decay };

// Parallel solver only: how the grid is split over MPI ranks, and how halos
// are exchanged (docs/design.md, "Domain decomposition" and "Halo exchange").
enum class Decomp { rows, two_d };
enum class Exchange { blocking, overlap };

// Every member has a default value (written after "="), so  Config c;  is
// already the configuration from config/defaults.json.
struct Config {
    // Grid and geometry
    int nx = 256;
    int ny = 256;
    double lx = 1.0;
    double ly = 1.0;
    int steps = 10000;

    // Physics and time stepping
    double kappa = 0.005;
    double cfl = 0.8;
    // std::optional holds either a double or nothing. Empty means
    // "choose dt automatically from cfl".
    std::optional<double> dt;

    // Double-gyre flow
    double flow_amplitude = 0.1;
    double flow_epsilon = 0.25;
    double flow_period = 10.0;

    // Initial condition
    InitialCase initial_case = InitialCase::pigment_drop;
    double background = 0.0;
    double amplitude = 1.0;
    double blob_x = 0.5;
    double blob_y = 0.5;
    double blob_width = 0.05;

    // Output: path of the final field (.csv, or .npy for large grids);
    // empty means no field output.
    std::string output;

    // Snapshots for watching the dynamics: every snapshot_interval steps
    // (0 = off) the field is written to snapshot_dir as .npy files, averaged
    // over snapshot_coarsen x snapshot_coarsen blocks of cells (1 = full size).
    int snapshot_interval = 0;
    std::string snapshot_dir = "outputs/snapshots";
    int snapshot_coarsen = 1;

    // Parallel execution (heat_parallel only). py, px: process grid rows and
    // columns for --decomp 2d; 0 means "choose automatically".
    Decomp decomposition = Decomp::rows;
    Exchange exchange = Exchange::blocking;
    int py = 0;
    int px = 0;
};

// The command-line help text.
std::string usage(const std::string& program);

// True if the arguments ask for help (-h or --help).
bool wants_help(int argc, char** argv);

// Read flags such as  --nx 128 --kappa 0.01  into a Config.
// Throws std::invalid_argument on an unknown flag or a missing value.
Config parse_args(int argc, char** argv);

// Names as written on the command line, e.g. "pigment_drop", "2d", "overlap".
std::string case_name(InitialCase c);
std::string decomp_name(Decomp d);
std::string exchange_name(Exchange e);

// True if any heat_parallel-only option differs from its default.
bool has_parallel_options(const Config& c);

// Throws std::invalid_argument if any parameter is outside its allowed range.
void validate(const Config& c);

// The flow parameters in the form velocity.hpp expects.
FlowParams flow_params(const Config& c);

// The time step: cfl / (a + d), or the explicit dt if given and stable.
double time_step(const Config& c);

// Peclet number U L / kappa, with U = max(u_max, v_max) and L = max(lx, ly).
double peclet(const Config& c);

// Grid Peclet number max(u_max dx, v_max dy) / kappa. Well above 2, upwind
// numerical diffusion is comparable to the physical diffusion.
double peclet_grid(const Config& c);
