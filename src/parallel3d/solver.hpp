// The MPI + OpenMP 3D pigment solver (model A of docs/design_3d.md): the 3D
// pigment in the 2D wind-driven flow of heat_parallel.

#pragma once

#include <mpi.h>

#include <vector>

#include "common/io.hpp"
#include "common3d/config3d.hpp"
#include "parallel/parallel_tools.hpp"  // configure_threads

struct Parallel3DResult {
    bool is_root = false;       // rank 0 of the solver's communicator
    RunSummary summary;         // complete on the root rank
    std::vector<double> field;  // final global nz x ny x nx field on the root rank if cfg.base.output is set
};

// Run the configured problem on all ranks of `world`, with the threads of
// configure_threads(cfg.base.threads) per rank. The ranks split x and y only
// (slabs: y; pencils: x and y); every rank holds all nz layers of its block.
// Collective: every rank must call it with the same config. Writes snapshots
// (from the root rank) if cfg.base.snapshot_interval > 0.
Parallel3DResult run_parallel_3d(const Config3D& c3, MPI_Comm world);
