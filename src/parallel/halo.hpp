// Halo exchange of the scalar field between neighbouring ranks.
// See docs/design.md, "Halo exchange".
//
// Rows are contiguous in memory and are sent directly from the field with an
// MPI vector datatype (halo rows of nx interior values). Columns are not
// contiguous: they are packed into buffers before sending and unpacked after
// receiving. The same scheme is used for blocking and overlapping exchange, so
// comparing the two changes nothing but the timing of the messages.
//
// Only interior cells are sent and only edge ghost cells are written; corner
// ghost cells are never read by the five-point update and are left alone.

#pragma once

#include <mpi.h>

#include <vector>

#include "common/grid.hpp"
#include "parallel/decomposition.hpp"

class HaloExchange {
public:
    HaloExchange(const ProcessLayout& layout, int ny, int nx);
    ~HaloExchange();
    HaloExchange(const HaloExchange&) = delete;
    HaloExchange& operator=(const HaloExchange&) = delete;

    // Blocking: MPI_Sendrecv in each direction. On return, every process
    // halo of `c` holds its neighbour's values.
    void exchange(Field& c);

    // Nonblocking: post all receives, pack the columns and post all sends.
    // Until finish(), do not modify the interior of `c` or read its process halos.
    void start(Field& c);

    // Wait for all messages of start() and unpack the received columns.
    void finish(Field& c);

private:
    void pack_columns(const Field& c);
    void unpack_columns(Field& c);

    MPI_Comm comm_;
    int north_;
    int south_;
    int west_;
    int east_;
    int ny_;
    int nx_;
    MPI_Datatype rows_;  // `halo` rows of nx interior values, stride = field stride
    std::vector<double> send_west_;
    std::vector<double> send_east_;
    std::vector<double> recv_west_;
    std::vector<double> recv_east_;
    MPI_Request requests_[8];
    int pending_ = 0;
};
