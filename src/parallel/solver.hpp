// The MPI + OpenMP solver. See src/README.md, "Parallel solver".

#pragma once

#include <mpi.h>

#include <vector>

#include "common/config.hpp"
#include "common/io.hpp"

struct ParallelResult {
    bool is_root = false;       // rank 0 of the solver's communicator
    RunSummary summary;         // complete on the root rank
    std::vector<double> field;  // final global ny x nx field on the root rank, empty elsewhere
};

// Run the configured problem on all ranks of `world`, with OMP_NUM_THREADS
// threads per rank. Collective: every rank must call it with the same config.
// Writes snapshots (from the root rank) if cfg.snapshot_interval > 0.
ParallelResult run_parallel(const Config& cfg, MPI_Comm world);
