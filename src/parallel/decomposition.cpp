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
    MPI_Comm_size(world, &size);

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
            MPI_Dims_create(size, 2, dims);  // as square as possible, dims[0] >= dims[1]
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
    MPI_Cart_create(world, 2, dims, periods, 1, &layout.comm);
    MPI_Comm_rank(layout.comm, &layout.rank);
    layout.size = size;
    layout.py = dims[0];
    layout.px = dims[1];

    int coords[2] = {0, 0};
    MPI_Cart_coords(layout.comm, layout.rank, 2, coords);
    layout.ry = coords[0];
    layout.rx = coords[1];

    // Dimension 0 runs down the rows: the source of a +1 shift is the rank above.
    MPI_Cart_shift(layout.comm, 0, 1, &layout.north, &layout.south);
    MPI_Cart_shift(layout.comm, 1, 1, &layout.west, &layout.east);

    Block& b = layout.block;
    split_evenly(c.ny, layout.py, layout.ry, b.ny, b.row_offset);
    split_evenly(c.nx, layout.px, layout.rx, b.nx, b.col_offset);
    b.wall_top = layout.ry == 0;
    b.wall_bottom = layout.ry == layout.py - 1;
    b.wall_left = layout.rx == 0;
    b.wall_right = layout.rx == layout.px - 1;
    return layout;
}
