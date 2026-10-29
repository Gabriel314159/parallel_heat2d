// Function bodies for parallel/decomposition.hpp.

#include "parallel/decomposition.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

void split_evenly(int n, int parts, int index, int& count, int& offset) {
    const int base = n / parts;
    const int extra = n % parts;
    count = base + (index < extra ? 1 : 0);
    offset = index * base + std::min(index, extra);
}

ProcessLayout make_layout(const Config& c, MPI_Comm world) {
    int size = 1;
    mpi_check(MPI_Comm_size(world, &size), "MPI_Comm_size");

    int dims[2] = {size, 1};  // {process rows, process columns}
    if (c.decomposition == Decomp::two_d) {
        if (c.py > 0) {
            if (c.py * c.px != size) {
                throw std::invalid_argument("--py " + std::to_string(c.py) + " x --px " + std::to_string(c.px) +
                                            " does not match " + std::to_string(size) + " MPI ranks");
            }
            dims[0] = c.py;
            dims[1] = c.px;
        } else {
            dims[0] = 0;
            dims[1] = 0;
            // As square as possible, dims[0] >= dims[1].
            mpi_check(MPI_Dims_create(size, 2, dims), "MPI_Dims_create");
            if (c.nx > c.ny) std::swap(dims[0], dims[1]);  // more process columns for a wide grid
        }
    }
    if (dims[0] > c.ny || dims[1] > c.nx) {
        throw std::invalid_argument("a " + std::to_string(dims[0]) + " x " + std::to_string(dims[1]) +
                                    " process grid leaves empty blocks on a " + std::to_string(c.ny) + " x " +
                                    std::to_string(c.nx) + " grid");
    }

    ProcessLayout layout;
    const int periods[2] = {0, 0};
    // reorder = 0: the ranks of the Cartesian communicator are the ranks of
    // `world`, so MPI cannot move a rank away from the core the launcher chose.
    mpi_check(MPI_Cart_create(world, 2, dims, periods, 0, &layout.comm), "MPI_Cart_create");
    mpi_check(MPI_Comm_set_errhandler(layout.comm, MPI_ERRORS_RETURN), "MPI_Comm_set_errhandler");
    mpi_check(MPI_Comm_rank(layout.comm, &layout.rank), "MPI_Comm_rank");
    layout.size = size;
    layout.py = dims[0];
    layout.px = dims[1];

    int coords[2] = {0, 0};
    mpi_check(MPI_Cart_coords(layout.comm, layout.rank, 2, coords), "MPI_Cart_coords");
    layout.ry = coords[0];
    layout.rx = coords[1];

    // Dimension 0 runs down the rows: the source of a +1 shift is the rank above.
    mpi_check(MPI_Cart_shift(layout.comm, 0, 1, &layout.north, &layout.south), "MPI_Cart_shift rows");
    mpi_check(MPI_Cart_shift(layout.comm, 1, 1, &layout.west, &layout.east), "MPI_Cart_shift columns");

    Block& b = layout.block;
    split_evenly(c.ny, layout.py, layout.ry, b.ny, b.row_offset);
    split_evenly(c.nx, layout.px, layout.rx, b.nx, b.col_offset);
    b.wall_top = layout.ry == 0;
    b.wall_bottom = layout.ry == layout.py - 1;
    b.wall_left = layout.rx == 0;
    b.wall_right = layout.rx == layout.px - 1;
    return layout;
}

int count_nodes(const ProcessLayout& layout) {
    // The ranks of one node form one shared-memory communicator. Counting its
    // rank 0 once per node gives the number of nodes.
    MPI_Comm shared = MPI_COMM_NULL;
    mpi_check(MPI_Comm_split_type(layout.comm, MPI_COMM_TYPE_SHARED, layout.rank, MPI_INFO_NULL, &shared),
              "shared-memory communicator");
    int shared_rank = 0;
    mpi_check(MPI_Comm_rank(shared, &shared_rank), "MPI_Comm_rank (shared-memory communicator)");
    const int leader = shared_rank == 0 ? 1 : 0;
    int nodes = 0;
    mpi_check(MPI_Allreduce(&leader, &nodes, 1, MPI_INT, MPI_SUM, layout.comm), "node count");
    MPI_Comm_free(&shared);
    return nodes;
}
