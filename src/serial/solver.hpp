// The serial reference solver. See src/README.md, "Serial solver".

#pragma once

#include <vector>

#include "common/config.hpp"
#include "common/io.hpp"

// The name the program reports (heat_serial_omp: see CMakeLists.txt).
#ifdef HEAT2D_SERIAL_OMP
inline constexpr const char* serial_solver_name = "heat_serial_omp";
#else
inline constexpr const char* serial_solver_name = "heat_serial";
#endif

struct SerialResult {
    RunSummary summary;
    std::vector<double> field;  // final interior field, ny x nx, row by row
};

// Run the configured problem on one core. Writes snapshots if
// cfg.snapshot_interval > 0; the final field is returned, not written.
SerialResult run_serial(const Config& cfg);
