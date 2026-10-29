// Output: final fields, snapshot series and the run summary.
//
// Formats
// - Final field (--output PATH): the pigment as CSV (ny rows of nx values, top
//   row first, 17 significant digits) or as .npy. Next to it, PATH with the
//   extension replaced by .json holds the run summary.
// - Snapshots (--snapshot-interval N): in the snapshot folder,
//     meta.json          the run configuration and dt,
//     index.csv          step,time,mass,energy,enstrophy,file,psi_file,vort_file
//                        (one line per snapshot; time is the flow time),
//     field_NNNNNN.npy   pigment, ny x nx cell values (float64),
//     psi_NNNNNN.npy     streamfunction on the (ny+1) x (nx+1) corners, walls included,
//     vort_NNNNNN.npy    vorticity on the same corners.
//   With --snapshot-coarsen K, the pigment is averaged over K x K blocks and
//   psi and the vorticity are sampled at every K-th corner.
// - Run summary: printed at the end of a run; the last line is
//     RESULT {...json...}
//   so scripts can find and parse it.

#pragma once

#include <fstream>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "common/config.hpp"
#include "common/grid.hpp"

struct RunSummary {
    std::string solver;  // "heat_serial" or "heat_parallel"
    int mpi_ranks = 1;
    int omp_threads = 1;      // threads per rank actually in use
    int nodes = 1;            // nodes the ranks run on
    int py = 1;
    int px = 1;
    std::string decomposition = "none";
    std::string exchange = "none";
    double dt = 0.0;
    int spinup_steps = 0;
    double t_release = 0.0;   // flow time when the pigment is released (end of spin-up)
    double t_end = 0.0;       // flow time at the end of the run
    double elapsed_s = 0.0;   // main time-stepping loop; snapshots excluded; max over ranks
    double spinup_s = 0.0;    // spin-up loop; max over ranks
    double poisson_s = 0.0;   // part of spin-up and main loop spent in the Poisson solve
    long poisson_iterations = 0;     // CG iterations over all flow steps
    int poisson_iterations_max = 0;  // the most in one solve
    long flow_steps = 0;             // flow steps (spin-up and main loop)
    double snapshot_s = 0.0;  // time spent gathering and writing snapshots
    double mass_initial = 0.0;
    double mass_final = 0.0;
    // |M_end - M_0| / (dx dy sum |c_0|): the relative mass change, also
    // meaningful when M_0 = 0 (cosine_decay).
    double mass_error = 0.0;
    double energy = 0.0;      // kinetic energy of the flow at the end
    double enstrophy = 0.0;   // enstrophy of the flow at the end
    double u_max_measured = 0.0;  // largest face speeds seen at the checks
    double v_max_measured = 0.0;
    double peclet = 0.0;      // with the estimated --max-speed
    double peclet_grid = 0.0; // with the largest measured speed
    // Only for cosine_decay; NaN otherwise (written as null).
    double cosine_discrete_error = std::numeric_limits<double>::quiet_NaN();
    double cosine_continuous_error = std::numeric_limits<double>::quiet_NaN();
    std::string compiler;
    std::string mpi_library = "none";
    std::string build_type;
};

// Compiler name and version, and the CMake build type this was built with.
std::string compiler_id();
std::string build_type();

// The configuration and summary as one line of JSON.
std::string summary_json(const Config& c, const RunSummary& s);

// A short human-readable report followed by the RESULT line.
void print_summary(std::ostream& out, const Config& c, const RunSummary& s);

// The interior cells of a field (no ghosts), ny x nx, row by row.
std::vector<double> copy_interior(const Field& field);

// A corner field of a block covering the whole grid as (ny+1) x (nx+1) values:
// corner rows 0..ny and columns 0..nx, walls included.
std::vector<double> copy_corners(const Field& field);

// The same, from the ny x nx entries (corner rows and columns 1..ny, 1..nx)
// gathered from all ranks: the wall corners of row 0 and column 0 are added as 0.
std::vector<double> corners_from_entries(const std::vector<double>& entries, int ny, int nx);

// Create the folder that will contain `path`, if needed.
void ensure_parent_directory(const std::string& path);

// `path` with its extension replaced, e.g. ("out/final.csv", ".json") -> "out/final.json".
std::string replace_extension(const std::string& path, const std::string& extension);

void write_field_csv(const std::string& path, const std::vector<double>& values, int ny, int nx);
void write_field_npy(const std::string& path, const std::vector<double>& values, int ny, int nx);
void write_text_file(const std::string& path, const std::string& text);

// The final field, as CSV or .npy depending on the extension of `path`.
void write_field(const std::string& path, const std::vector<double>& values, int ny, int nx);

// Average an ny x nx field over k x k blocks (k must divide ny and nx). Block
// averages keep the total mass: sum(coarse) * k^2 = sum(fine).
std::vector<double> coarsen_field(const std::vector<double>& values, int ny, int nx, int k);

// Every k-th value of an (ny+1) x (nx+1) corner field (k must divide ny and nx).
std::vector<double> subsample_corners(const std::vector<double>& values, int ny, int nx, int k);

// One snapshot, all fields on the whole grid.
struct Snapshot {
    int step = 0;           // steps since the pigment was released
    double time = 0.0;      // flow time
    double mass = 0.0;
    double energy = 0.0;
    double enstrophy = 0.0;
    std::vector<double> pigment;  // ny x nx
    std::vector<double> psi;      // (ny+1) x (nx+1)
    std::vector<double> zeta;     // (ny+1) x (nx+1)
};

// Writes the snapshot folder described at the top of this file.
class SnapshotWriter {
public:
    SnapshotWriter(const Config& c, double dt, double t_release, const std::string& solver);
    void write(const Snapshot& s);

private:
    std::string dir_;
    int ny_;
    int nx_;
    int coarsen_;
    std::ofstream index_;
};
