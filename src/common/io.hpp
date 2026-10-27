// Output: final fields, snapshot series and the run summary.
//
// Formats
// - Final field (--output PATH): CSV, ny rows of nx values, top row first,
//   17 significant digits. Next to it, PATH with the extension replaced by
//   .json holds the run summary.
// - Snapshots (--snapshot-interval N): in the snapshot folder,
//     meta.json          the run configuration and dt,
//     index.csv          step,time,mass,file  (one line per snapshot),
//     field_000000.npy   ny x nx float64 arrays (numpy.load reads them).
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
    int omp_threads = 1;
    int py = 1;
    int px = 1;
    std::string decomposition = "none";
    std::string exchange = "none";
    double dt = 0.0;
    double t_end = 0.0;
    double elapsed_s = 0.0;   // time-stepping loop; snapshots excluded; max over ranks
    double snapshot_s = 0.0;  // time spent gathering and writing snapshots
    double mass_initial = 0.0;
    double mass_final = 0.0;
    // |M_end - M_0| / (dx dy sum |c_0|): the relative mass change, also
    // meaningful when M_0 = 0 (cosine_decay).
    double mass_error = 0.0;
    double peclet = 0.0;
    double peclet_grid = 0.0;
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

// Writes the snapshot folder described at the top of this file. `values` is
// the full ny x nx field; with --snapshot-coarsen K it is stored as block
// averages, (ny/K) x (nx/K). `mass` is always that of the full field.
class SnapshotWriter {
public:
    SnapshotWriter(const Config& c, double dt, const std::string& solver);
    void write(int step, const std::vector<double>& values, double mass);

private:
    std::string dir_;
    double dt_;
    int ny_;
    int nx_;
    int coarsen_;
    std::ofstream index_;
};
