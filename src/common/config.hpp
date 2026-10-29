// Run configuration: command-line parsing, validation and the global time step.
// Flags and defaults follow config/README.md and config/defaults.json.

#pragma once

#include <optional>
#include <string>

// "enum class" is a type with a fixed set of named values. A variable of this
// type can only hold one of them, so a typo like "pigment_dorp" cannot sneak
// through as it could with a plain string.
enum class InitialCase { pigment_drop, constant, cosine_decay };

// The two forms of wind forcing (docs/design.md, "The wind-driven flow").
enum class Wind { steady, periodic };

// Parallel solver only: how the grid is split over MPI ranks, and how halos
// are exchanged (docs/design.md, "Domain decomposition" and "Halo exchange").
enum class Decomp { rows, two_d };
enum class Exchange { blocking, overlap };

// Every member has a default value (written after "="), so  Config c;  is
// already the configuration from config/defaults.json.
struct Config {
    // Grid and geometry: a 1 x 2 basin with square cells.
    int nx = 256;
    int ny = 512;
    double lx = 1.0;
    double ly = 2.0;
    int steps = 50000;

    // Pigment and time stepping
    double kappa = 1e-4;
    double cfl = 0.8;
    // std::optional holds either a double or nothing. Empty means
    // "choose dt automatically from cfl and max_speed".
    std::optional<double> dt;

    // Wind-driven flow (docs/design.md, "The wind-driven flow"). With
    // nu = 1e-3 and F0 = 0.1 the flow reaches speeds of about 1 (Re about
    // 1000) and settles into two gyres over a spin-up of about 100.
    double viscosity = 1e-3;
    Wind wind = Wind::steady;
    double wind_amplitude = 0.1;
    double wind_alpha = 0.0;
    double wind_period = 1.0;
    double spinup = 100.0;  // flow-only time before the pigment is released

    // The fixed time step is derived from this estimate of the largest flow
    // speed, and every check_interval steps (0 = never) the actual speeds are
    // checked against the stability limits.
    double max_speed = 2.0;
    int check_interval = 100;

    // Poisson solver (conjugate gradients): stop when the residual is below
    // poisson_tol times the norm of the right-hand side; give up after
    // poisson_max_iter iterations.
    double poisson_tol = 1e-10;
    int poisson_max_iter = 20000;

    // Initial condition
    InitialCase initial_case = InitialCase::pigment_drop;
    double background = 0.0;
    double amplitude = 1.0;
    double blob_x = 0.2;
    double blob_y = 1.0;
    double blob_width = 0.05;

    // Output: path of the final field (.csv, or .npy for large grids);
    // empty means no field output.
    std::string output;

    // Snapshots for watching the dynamics: every snapshot_interval steps
    // (0 = off) the fields are written to snapshot_dir as .npy files,
    // coarsened by snapshot_coarsen (1 = full size).
    int snapshot_interval = 0;
    std::string snapshot_dir = "outputs/snapshots";
    int snapshot_coarsen = 1;

    // Parallel execution (heat_parallel only). py, px: process grid rows and
    // columns for --decomp 2d; 0 means "choose automatically".
    Decomp decomposition = Decomp::rows;
    Exchange exchange = Exchange::blocking;
    int py = 0;
    int px = 0;

    // OpenMP threads per rank. Empty means "use OMP_NUM_THREADS, or one
    // thread if that is not set". heat_serial accepts 1 only.
    std::optional<int> threads;
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
std::string wind_name(Wind w);
std::string decomp_name(Decomp d);
std::string exchange_name(Exchange e);

// True if any heat_parallel-only option differs from its default.
bool has_parallel_options(const Config& c);

// Throws std::invalid_argument if any parameter is outside its allowed range.
void validate(const Config& c);

// The time step: cfl times the largest stable step for speed max_speed, or the
// explicit dt if given and stable for that speed (see common/stability.hpp).
double time_step(const Config& c);

// Number of flow-only steps before the pigment is released: spinup / dt, rounded.
int spinup_steps(const Config& c, double dt);

// Peclet number U L / kappa, with U = max_speed and L = max(lx, ly).
double peclet(const Config& c);

// Grid Peclet number U max(dx, dy) / kappa for a speed U. Well above 2,
// upwind numerical diffusion is comparable to the physical diffusion.
double peclet_grid(const Config& c, double speed);
