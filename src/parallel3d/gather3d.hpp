// Collect pieces of a global 3D array on rank 0, for output only (never in
// the timed loop): the whole 3D pigment, its surface layer, a vertical
// section, and the 2D psi and vorticity of the snapshots.

#pragma once

#include <mpi.h>

#include <climits>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "parallel/decomposition.hpp"

// The part of a global nk x ni x nj array that one rank holds: offsets and
// sizes along (k, i, j). An empty box (a size of 0) holds nothing.
struct Box {
    int k0 = 0;
    int i0 = 0;
    int j0 = 0;
    int nk = 0;
    int ni = 0;
    int nj = 0;
    std::size_t count() const {
        return static_cast<std::size_t>(nk) * static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
    }
};

// Collective. Every rank passes its box and its values (k slowest, j
// fastest); rank 0 returns the global array, the other ranks an empty vector.
// Positions no box covers are 0. MPI counts are int, so arrays of 2^31 values
// or more are refused; write snapshots (one layer and one section) instead.
inline std::vector<double> gather_boxes(MPI_Comm comm, const Box& mine, const std::vector<double>& values, int nk,
                                        int ni, int nj) {
    const long long total = static_cast<long long>(nk) * ni * nj;
    if (total >= INT_MAX) {
        throw std::runtime_error("cannot gather " + std::to_string(total) +
                                 " values on one rank (MPI counts are int); write snapshots instead of --output");
    }
    if (values.size() != mine.count()) throw std::logic_error("gather_boxes: values do not match the box");
    int rank = 0;
    int size = 1;
    mpi_check(MPI_Comm_rank(comm, &rank), "MPI_Comm_rank (gather)");
    mpi_check(MPI_Comm_size(comm, &size), "MPI_Comm_size (gather)");

    const int box[6] = {mine.k0, mine.i0, mine.j0, mine.nk, mine.ni, mine.nj};
    std::vector<int> boxes(rank == 0 ? static_cast<std::size_t>(6) * size : 0);
    mpi_check(MPI_Gather(box, 6, MPI_INT, rank == 0 ? boxes.data() : nullptr, 6, MPI_INT, 0, comm), "gather boxes");

    std::vector<int> counts(rank == 0 ? size : 0);
    std::vector<int> displacements(rank == 0 ? size : 0);
    int offset = 0;
    if (rank == 0) {
        for (int r = 0; r < size; ++r) {
            const int* b = &boxes[static_cast<std::size_t>(6) * r];
            counts[r] = b[3] * b[4] * b[5];
            displacements[r] = offset;
            offset += counts[r];
        }
    }
    std::vector<double> received(rank == 0 ? static_cast<std::size_t>(offset) : 0);
    mpi_check(MPI_Gatherv(values.data(), static_cast<int>(values.size()), MPI_DOUBLE, received.data(),
                          rank == 0 ? counts.data() : nullptr, rank == 0 ? displacements.data() : nullptr, MPI_DOUBLE,
                          0, comm),
              "gather values");
    if (rank != 0) return {};

    std::vector<double> global(static_cast<std::size_t>(total), 0.0);
    for (int r = 0; r < size; ++r) {
        const int* b = &boxes[static_cast<std::size_t>(6) * r];
        const double* piece = received.data() + displacements[r];
        std::size_t at = 0;
        for (int k = 0; k < b[3]; ++k) {
            for (int i = 0; i < b[4]; ++i) {
                for (int j = 0; j < b[5]; ++j) {
                    const std::size_t to = (static_cast<std::size_t>(b[0] + k) * ni + static_cast<std::size_t>(b[1] + i)) *
                                               static_cast<std::size_t>(nj) +
                                           static_cast<std::size_t>(b[2] + j);
                    global[to] = piece[at++];
                }
            }
        }
    }
    return global;
}
