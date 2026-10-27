// Domain decomposition: which block of the global grid each MPI rank owns,
// and who its neighbours are. See docs/design.md, "Domain decomposition".

#pragma once

#include <mpi.h>

#include "common/config.hpp"
#include "common/grid.hpp"

// Split n items into `parts` nearly equal parts: the first n % parts parts get
// one extra item. Returns the size and the offset of part `index`.
void split_evenly(int n, int parts, int index, int& count, int& offset);

struct ProcessLayout {
    // Cartesian communicator. Use its rank numbers everywhere: with rank
    // reordering they may differ from MPI_COMM_WORLD.
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0;
    int size = 1;
    int py = 1;  // process grid rows
    int px = 1;  // process grid columns
    int ry = 0;  // this rank's process row, 0 at the top
    int rx = 0;  // this rank's process column, 0 at the left
    // Neighbouring ranks, or MPI_PROC_NULL on a physical wall.
    int north = MPI_PROC_NULL;
    int south = MPI_PROC_NULL;
    int west = MPI_PROC_NULL;
    int east = MPI_PROC_NULL;
    Block block;  // this rank's cells
};

// Build the process grid for --decomp rows (p x 1) or 2d (--py x --px, or
// chosen with MPI_Dims_create). Throws std::invalid_argument if the grid does
// not fit, e.g. more process rows than cell rows. Collective over `world`.
ProcessLayout make_layout(const Config& c, MPI_Comm world);

// Frees the communicator of a layout when it goes out of scope.
class LayoutGuard {
public:
    explicit LayoutGuard(ProcessLayout& layout) : layout_(layout) {}
    ~LayoutGuard() {
        if (layout_.comm != MPI_COMM_NULL) MPI_Comm_free(&layout_.comm);
    }
    LayoutGuard(const LayoutGuard&) = delete;
    LayoutGuard& operator=(const LayoutGuard&) = delete;

private:
    ProcessLayout& layout_;
};
