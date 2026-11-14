// Halo exchange of the 3D pigment between neighbouring ranks
// (docs/design_3d.md, "Model A", "Halo exchange").
//
// Model A splits only x and y (slabs or pencils), so every rank owns whole
// columns of water and has at most four neighbours: north, south, west and
// east. Each face it exchanges is a 2D slab, one row or one column in every
// layer. The surface and the bottom are walls of every rank, filled locally.
//
// - North/south faces: row 1 (or ny) of every layer is nx contiguous values
//   per layer, one layer apart, so it is sent directly from the field with an
//   MPI vector datatype (nz blocks of nx values, stride = one stored layer).
// - West/east faces: column 1 (or nx) of every row of every layer is not
//   contiguous; it is packed into a buffer of nz x ny values, as the 2D
//   columns are.
//
// As in 2D, only interior cells are sent and only face ghosts are written; the
// seven-point update reads no edge or corner ghosts.

#pragma once

#include <mpi.h>

#include <vector>

#include "common/grid.hpp"
#include "common3d/field3d.hpp"
#include "parallel/decomposition.hpp"

class Halo3D {
public:
    Halo3D(const ProcessLayout& layout, int nz, int ny, int nx);
    ~Halo3D();
    Halo3D(const Halo3D&) = delete;
    Halo3D& operator=(const Halo3D&) = delete;

    // Blocking: MPI_Sendrecv in each direction. On return, every process
    // halo of `c` holds its neighbour's values.
    void exchange(Field3D& c);

    // Nonblocking: post all receives, pack the columns and post all sends.
    // Until finish(), do not modify the interior of `c` or read its process halos.
    void start(Field3D& c);

    // Wait for all messages of start() and unpack the received columns.
    void finish(Field3D& c);

private:
    void pack_columns(const Field3D& c);
    void unpack_columns(Field3D& c);

    MPI_Comm comm_;
    int north_;
    int south_;
    int west_;
    int east_;
    int nz_;
    int ny_;
    int nx_;
    MPI_Datatype rows_;  // one row of nx interior values in every layer
    std::vector<double> send_west_;
    std::vector<double> send_east_;
    std::vector<double> recv_west_;
    std::vector<double> recv_east_;
    MPI_Request requests_[8];
    int pending_ = 0;
};
