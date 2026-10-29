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

// Set the number of OpenMP threads of this rank: `requested` (--threads) if it
// is positive, else OMP_NUM_THREADS if that is set, else 1. Dynamic adjustment
// is switched off, so the number is what was asked for. Returns the number of
// threads in use. Throws if more than one thread is requested in a build
// without OpenMP. Call it before the first OpenMP region.
int configure_threads(int requested);

// Run the configured problem on all ranks of `world`, with the threads of
// configure_threads(cfg.threads) per rank. Collective: every rank must call it
// with the same config. Writes snapshots (from the root rank) if
// cfg.snapshot_interval > 0.
ParallelResult run_parallel(const Config& cfg, MPI_Comm world);
