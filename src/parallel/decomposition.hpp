// Domain decomposition: which block of the global grid each MPI rank owns,
// and who its neighbours are. See docs/design.md, "Domain decomposition".

#pragma once

#include <mpi.h>

#include <cstddef>
#include <stdexcept>
#include <string>

#include "common/config.hpp"
#include "common/grid.hpp"

// Turn the return code of an MPI call into an exception that names the call.
// The communicators use MPI_ERRORS_RETURN, so a failing call returns a code
// instead of aborting on the spot. The exception ends in main(), which prints
// the rank and stops all ranks with MPI_Abort. Wrap every MPI call in this,
// except in destructors, which must not throw.
inline void mpi_check(int code, const char* operation) {
    if (code != MPI_SUCCESS) {
        char text[MPI_MAX_ERROR_STRING];
        int length = 0;
        MPI_Error_string(code, text, &length);
        throw std::runtime_error(std::string(operation) + ": " +
                                 std::string(text, static_cast<std::size_t>(length)));
    }
}

// Split n items into `parts` nearly equal parts: the first n % parts parts get
// one extra item. Returns the size and the offset of part `index`.
void split_evenly(int n, int parts, int index, int& count, int& offset);

struct ProcessLayout {
    // Cartesian communicator, created without rank reordering: its rank
    // numbers are those of the communicator it was made from (MPI_COMM_WORLD),
    // so the placement of ranks on cores stays as the launcher set it. Used
    // for all communication after the layout is made. Its error handler is
    // MPI_ERRORS_RETURN (see mpi_check).
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

// The number of nodes the ranks of a layout run on: ranks that can share
// memory (MPI_COMM_TYPE_SHARED) form one node. Collective.
int count_nodes(const ProcessLayout& layout);

// Frees the communicator of a layout when it goes out of scope.
class LayoutGuard {
public:
    explicit LayoutGuard(ProcessLayout& layout) : layout_(layout) {}
    ~LayoutGuard() {
        // No mpi_check: a destructor must not throw.
        if (layout_.comm != MPI_COMM_NULL) MPI_Comm_free(&layout_.comm);
    }
    LayoutGuard(const LayoutGuard&) = delete;
    LayoutGuard& operator=(const LayoutGuard&) = delete;

private:
    ProcessLayout& layout_;
};
