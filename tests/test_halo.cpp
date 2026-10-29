// MPI test for src/parallel/: the blocking, nonblocking and corner halo
// exchanges, and the setup of the process layout. Every cell holds a value
// encoding its global position; after an exchange, every ghost cell next to
// another rank must hold the value of the cell it mirrors. The layout must
// keep the world ranks (no reordering), return MPI errors instead of aborting,
// and count the nodes. Run by CTest on several rank counts.

#include <mpi.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include "common/config.hpp"
#include "common/corners.hpp"
#include "common/grid.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"

namespace {

double code(int ig, int jg) {
    return 1000.0 * ig + jg;
}

// Fill the interior with codes of global positions and the ghosts with -1.
void fill(Field& f, const Block& b) {
    f.fill(-1.0);
    for (int i = 1; i <= b.ny; ++i) {
        for (int j = 1; j <= b.nx; ++j) f(i, j) = code(b.row_offset + i, b.col_offset + j);
    }
}

// Count wrong ghost cells: edges always; with corners, also the four corner ghosts.
int wrong_ghosts(const Field& f, const ProcessLayout& l, bool corners) {
    const Block& b = l.block;
    int wrong = 0;
    auto expect = [&](int i, int j) {
        if (f(i, j) != code(b.row_offset + i, b.col_offset + j)) ++wrong;
    };
    const bool n = l.north != MPI_PROC_NULL, s = l.south != MPI_PROC_NULL;
    const bool w = l.west != MPI_PROC_NULL, e = l.east != MPI_PROC_NULL;
    for (int j = 1; j <= b.nx; ++j) {
        if (n) expect(0, j);
        if (s) expect(b.ny + 1, j);
    }
    for (int i = 1; i <= b.ny; ++i) {
        if (w) expect(i, 0);
        if (e) expect(i, b.nx + 1);
    }
    if (corners) {
        if (n && w) expect(0, 0);
        if (n && e) expect(0, b.nx + 1);
        if (s && w) expect(b.ny + 1, 0);
        if (s && e) expect(b.ny + 1, b.nx + 1);
    }
    return wrong;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int failures = 0;

    // mpi_check: silent on success; on an error code, throws and names the operation.
    bool silent_on_success = true;
    try {
        mpi_check(MPI_SUCCESS, "demo operation");
    } catch (...) {
        silent_on_success = false;
    }
    bool names_operation = false;
    try {
        mpi_check(MPI_ERR_OTHER, "demo operation");
    } catch (const std::runtime_error& e) {
        names_operation = std::string(e.what()).rfind("demo operation: ", 0) == 0;
    }
    if (rank == 0) {
        std::cout << (silent_on_success ? "  ok    " : "  FAIL  ") << "mpi_check does nothing on MPI_SUCCESS\n";
        std::cout << (names_operation ? "  ok    " : "  FAIL  ")
                  << "mpi_check throws on an error and names the operation\n";
        failures += !silent_on_success + !names_operation;
    }

    for (Decomp d : {Decomp::rows, Decomp::two_d}) {
        Config c;
        c.nx = 23;
        c.ny = 17;
        c.decomposition = d;
        ProcessLayout layout = make_layout(c, MPI_COMM_WORLD);
        const LayoutGuard guard(layout);
        HaloExchange halos(layout, layout.block.ny, layout.block.nx);
        Field f(layout.block.ny, layout.block.nx);

        // Layout setup, counted per rank: the Cartesian rank is the world rank
        // (no reordering); errors are returned, not fatal; the node count is
        // between 1 and the number of ranks, and 1 for a single rank.
        int setup[3] = {0, 0, 0};
        setup[0] = layout.rank != rank;
        MPI_Errhandler handler = MPI_ERRHANDLER_NULL;
        MPI_Comm_get_errhandler(layout.comm, &handler);
        setup[1] = handler != MPI_ERRORS_RETURN;
        MPI_Errhandler_free(&handler);
        const int nodes = count_nodes(layout);
        setup[2] = nodes < 1 || nodes > layout.size || (layout.size == 1 && nodes != 1);

        int wrong[3] = {0, 0, 0};
        fill(f, layout.block);
        halos.exchange(f);
        wrong[0] = wrong_ghosts(f, layout, false);
        fill(f, layout.block);
        halos.start(f);
        halos.finish(f);
        wrong[1] = wrong_ghosts(f, layout, false);
        fill(f, layout.block);
        halos.exchange_with_corners(f);
        wrong[2] = wrong_ghosts(f, layout, true);

        int total[3] = {0, 0, 0};
        MPI_Reduce(wrong, total, 3, MPI_INT, MPI_SUM, 0, layout.comm);
        int setup_total[3] = {0, 0, 0};
        MPI_Reduce(setup, setup_total, 3, MPI_INT, MPI_SUM, 0, layout.comm);
        if (rank == 0) {
            const std::string name = decomp_name(d) + " on " + std::to_string(layout.py) + " x " + std::to_string(layout.px);
            const char* setup_labels[3] = {"ranks equal the world ranks (no reordering)",
                                           "the Cartesian communicator returns MPI errors",
                                           "the node count is between 1 and the number of ranks"};
            for (int k = 0; k < 3; ++k) {
                std::cout << (setup_total[k] == 0 ? "  ok    " : "  FAIL  ") << name << ": " << setup_labels[k] << " ("
                          << setup_total[k] << " ranks wrong), nodes = " << nodes << "\n";
                failures += setup_total[k] != 0;
            }
            const char* labels[3] = {"blocking exchange", "nonblocking exchange", "exchange with corners"};
            for (int k = 0; k < 3; ++k) {
                std::cout << (total[k] == 0 ? "  ok    " : "  FAIL  ") << name << ": " << labels[k] << " ("
                          << total[k] << " wrong ghost cells)\n";
                failures += total[k] != 0;
            }
        }
    }

    if (rank == 0) std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    MPI_Bcast(&failures, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return failures == 0 ? 0 : 1;
}
