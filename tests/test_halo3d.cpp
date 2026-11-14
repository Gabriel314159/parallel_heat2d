// MPI test for src/parallel3d/halo3d.cpp: the halo exchange of the 3D pigment,
// for slabs and pencils, blocking and nonblocking. Every cell holds a value
// encoding its global position (layer, row, column); after an exchange, every
// face ghost next to another rank must hold the value of the cell it mirrors,
// in every layer, and no ghost on a wall (the side walls, the surface and the
// bottom) may be touched. Run by CTest on several rank counts.

#include <mpi.h>

#include <iostream>
#include <string>

#include "common/config.hpp"
#include "common/grid.hpp"
#include "common3d/config3d.hpp"
#include "common3d/field3d.hpp"
#include "parallel/decomposition.hpp"
#include "parallel3d/halo3d.hpp"

namespace {

double code(int k, int ig, int jg) {
    return 1e6 * k + 1000.0 * ig + jg;
}

// The interior holds the codes of global positions, every ghost -1.
void fill(Field3D& f, const Block& b) {
    f.fill(-1.0);
    for (int k = 1; k <= f.nz(); ++k) {
        for (int i = 1; i <= b.ny; ++i) {
            for (int j = 1; j <= b.nx; ++j) f(k, i, j) = code(k, b.row_offset + i, b.col_offset + j);
        }
    }
}

// Ghosts next to another rank must hold the neighbour's value; all other
// ghosts must still be -1 (walls are filled by the solver, not the exchange).
int wrong_ghosts(const Field3D& f, const ProcessLayout& l) {
    const Block& b = l.block;
    const bool n = l.north != MPI_PROC_NULL, s = l.south != MPI_PROC_NULL;
    const bool w = l.west != MPI_PROC_NULL, e = l.east != MPI_PROC_NULL;
    int wrong = 0;
    auto expect = [&](int k, int i, int j, bool neighbour) {
        const double want = neighbour ? code(k, b.row_offset + i, b.col_offset + j) : -1.0;
        if (f(k, i, j) != want) ++wrong;
    };
    for (int k = 1; k <= f.nz(); ++k) {
        for (int j = 1; j <= b.nx; ++j) {
            expect(k, 0, j, n);
            expect(k, b.ny + 1, j, s);
        }
        for (int i = 1; i <= b.ny; ++i) {
            expect(k, i, 0, w);
            expect(k, i, b.nx + 1, e);
        }
    }
    // The surface and bottom ghost layers belong to no neighbour.
    for (int i = 0; i <= b.ny + 1; ++i) {
        for (int j = 0; j <= b.nx + 1; ++j) {
            if (f(0, i, j) != -1.0 || f(f.nz() + 1, i, j) != -1.0) ++wrong;
        }
    }
    return wrong;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int failures = 0;

    for (Decomp d : {Decomp::rows, Decomp::two_d}) {
        for (int nz : {1, 4}) {
            Config c;
            c.nx = 23;
            c.ny = 17;
            c.decomposition = d;
            ProcessLayout layout = make_layout(c, MPI_COMM_WORLD);
            const LayoutGuard guard(layout);
            Halo3D halos(layout, nz, layout.block.ny, layout.block.nx);
            Field3D f(nz, layout.block.ny, layout.block.nx);

            int wrong[2] = {0, 0};
            fill(f, layout.block);
            halos.exchange(f);
            wrong[0] = wrong_ghosts(f, layout);
            fill(f, layout.block);
            halos.start(f);
            halos.finish(f);
            wrong[1] = wrong_ghosts(f, layout);
            // Twice in a row, as in the time loop.
            halos.start(f);
            halos.finish(f);
            wrong[1] += wrong_ghosts(f, layout);

            int total[2] = {0, 0};
            MPI_Reduce(wrong, total, 2, MPI_INT, MPI_SUM, 0, layout.comm);
            if (rank == 0) {
                const std::string name = decomp_name_3d(d) + " on " + std::to_string(layout.py) + " x " +
                                         std::to_string(layout.px) + ", " + std::to_string(nz) + " layers";
                const char* labels[2] = {"blocking exchange", "nonblocking exchange"};
                for (int k = 0; k < 2; ++k) {
                    std::cout << (total[k] == 0 ? "  ok    " : "  FAIL  ") << name << ": " << labels[k] << " ("
                              << total[k] << " wrong ghost cells)\n";
                    failures += total[k] != 0;
                }
            }
        }
    }

    if (rank == 0) std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    MPI_Bcast(&failures, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return failures == 0 ? 0 : 1;
}
