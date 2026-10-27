// The serial reference solver. See src/README.md, "Serial solver".

#pragma once

#include <vector>

#include "common/config.hpp"
#include "common/io.hpp"

struct SerialResult {
    RunSummary summary;
    std::vector<double> field;  // final interior field, ny x nx, row by row
};

// Run the configured problem on one core. Writes snapshots if
// cfg.snapshot_interval > 0; the final field is returned, not written.
SerialResult run_serial(const Config& cfg);
