// MPI test for src/parallel/: the blocking, nonblocking and corner halo
// exchanges, and the setup of the process layout. Every cell holds a value
// encoding its global position; after an exchange, every ghost cell next to
// another rank must hold the value of the cell it mirrors. The layout must
// keep the world ranks (no reordering), return MPI errors instead of aborting,
// and count the nodes. The CG Poisson solver, with blocking and with
// overlapping exchange, must give the serial solution. Run by CTest on several
// rank counts.

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "common/config.hpp"
#include "common/corners.hpp"
#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/poisson.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"
#include "parallel/poisson_comm.hpp"

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

// The right-hand side of the Poisson test at the corner with global indices (ig, jg).
double rhs(int ig, int jg) {
    return std::sin(0.37 * ig + 0.11) * std::cos(0.23 * jg - 0.4) + 0.001 * ig * jg;
}

// The overlapping exchange, but the process halos are set to NaN before the
// messages are started: a halo value read before finish() makes the result NaN.
class PoisonedHalos : public PoissonComm {
public:
    PoisonedHalos(const MpiComm& comm, const ProcessLayout& layout) : comm_(comm), layout_(layout) {}
    void sum(double* values, int count) const override { comm_.sum(values, count); }
    void exchange(Field& f) const override { comm_.exchange(f); }
    bool overlap() const override { return true; }
    void start(Field& f) const override {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const Block& b = layout_.block;
        for (int j = 1; j <= b.nx; ++j) {
            if (layout_.north != MPI_PROC_NULL) f(0, j) = nan;
            if (layout_.south != MPI_PROC_NULL) f(b.ny + 1, j) = nan;
        }
        for (int i = 1; i <= b.ny; ++i) {
            if (layout_.west != MPI_PROC_NULL) f(i, 0) = nan;
            if (layout_.east != MPI_PROC_NULL) f(i, b.nx + 1) = nan;
        }
        comm_.start(f);
    }
    void finish(Field& f) const override { comm_.finish(f); }

private:
    const MpiComm& comm_;
    const ProcessLayout& layout_;
};

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

    // The CG Poisson solver, blocking and overlapping (with poisoned halos),
    // against the serial solution of the whole grid, which every rank
    // computes. On the 6 x 6 grid some blocks have no inner corners.
    const int grids[2][2] = {{23, 17}, {6, 6}};
    for (const auto& grid : grids) {
        for (Decomp d : {Decomp::rows, Decomp::two_d}) {
            Config c;
            c.nx = grid[0];
            c.ny = grid[1];
            c.decomposition = d;
            const double dx = c.lx / c.nx;
            const double dy = c.ly / c.ny;

            const Block whole = whole_domain(c);
            const CornerRange all = unknown_corners(whole);
            Field zeta_all(c.ny, c.nx), psi_all(c.ny, c.nx);
            for (int i = all.i_first; i <= all.i_last; ++i) {
                for (int j = all.j_first; j <= all.j_last; ++j) zeta_all(i, j) = rhs(i, j);
            }
            PoissonCG serial(whole, dx, dy, 1e-12, 100000);
            serial.solve(zeta_all, psi_all, SerialComm());
            double scale = 0.0;
            for (int i = all.i_first; i <= all.i_last; ++i) {
                for (int j = all.j_first; j <= all.j_last; ++j) scale = std::max(scale, std::abs(psi_all(i, j)));
            }

            ProcessLayout layout = make_layout(c, MPI_COMM_WORLD);
            const LayoutGuard guard(layout);
            const Block& b = layout.block;
            const CornerRange u = unknown_corners(b);
            HaloExchange halos(layout, b.ny, b.nx);
            Field zeta(b.ny, b.nx), psi_blocking(b.ny, b.nx), psi_overlap(b.ny, b.nx);
            for (int i = u.i_first; i <= u.i_last; ++i) {
                for (int j = u.j_first; j <= u.j_last; ++j) zeta(i, j) = rhs(b.row_offset + i, b.col_offset + j);
            }
            MpiComm blocking_comm(layout.comm, halos, false);
            MpiComm overlap_comm(layout.comm, halos, true);
            const PoisonedHalos poisoned(overlap_comm, layout);
            // A solve that does not converge throws on all ranks together
            // (they see the same global sums); count it as wrong.
            CgStats stats[2];
            auto solve = [&](Field& psi, const PoissonComm& comm, CgStats& result) {
                PoissonCG solver(b, dx, dy, 1e-12, 5000);
                try {
                    result = solver.solve(zeta, psi, comm);
                    return 0;
                } catch (const std::runtime_error&) {
                    return 1;
                }
            };
            const int failed[2] = {solve(psi_blocking, blocking_comm, stats[0]), solve(psi_overlap, poisoned, stats[1])};

            // Differences from the serial psi, relative to its largest value;
            // a NaN counts as wrong.
            int wrong[2] = {0, 0};
            double largest[2] = {0.0, 0.0};
            for (int i = u.i_first; i <= u.i_last; ++i) {
                for (int j = u.j_first; j <= u.j_last; ++j) {
                    const double reference = psi_all(b.row_offset + i, b.col_offset + j);
                    const double diff[2] = {std::abs(psi_blocking(i, j) - reference) / scale,
                                            std::abs(psi_overlap(i, j) - reference) / scale};
                    for (int k = 0; k < 2; ++k) {
                        wrong[k] += !(diff[k] <= 1e-9);
                        largest[k] = std::max(largest[k], diff[k]);
                    }
                }
            }
            // With more than one rank, the timers must have recorded time.
            const CommTimes tb = blocking_comm.times();
            const CommTimes to = overlap_comm.times();
            int untimed = layout.size > 1 && !(tb.halo > 0 && tb.sum > 0 && to.halo > 0 && to.sum > 0);

            int wrong_total[2] = {0, 0};
            double largest_all[2] = {0.0, 0.0};
            int untimed_total = 0;
            MPI_Reduce(wrong, wrong_total, 2, MPI_INT, MPI_SUM, 0, layout.comm);
            MPI_Reduce(largest, largest_all, 2, MPI_DOUBLE, MPI_MAX, 0, layout.comm);
            MPI_Reduce(&untimed, &untimed_total, 1, MPI_INT, MPI_SUM, 0, layout.comm);
            if (rank == 0) {
                const std::string name = decomp_name(d) + " on " + std::to_string(layout.py) + " x " +
                                         std::to_string(layout.px) + ", " + std::to_string(c.nx) + " x " +
                                         std::to_string(c.ny) + " cells";
                const char* labels[2] = {"CG with blocking exchange", "CG with overlapping exchange, poisoned halos"};
                for (int k = 0; k < 2; ++k) {
                    const bool ok = failed[k] == 0 && wrong_total[k] == 0;
                    std::cout << (ok ? "  ok    " : "  FAIL  ") << name << ": " << labels[k] << " gives the serial psi (";
                    if (failed[k]) {
                        std::cout << "the solve did not converge)\n";
                    } else {
                        std::cout << stats[k].iterations << " iterations, largest relative difference "
                                  << largest_all[k] << ", " << wrong_total[k] << " corners wrong)\n";
                    }
                    failures += !ok;
                }
                std::cout << (untimed_total == 0 ? "  ok    " : "  FAIL  ") << name
                          << ": the CG's halo exchanges and global sums are timed\n";
                failures += untimed_total != 0;
            }
        }
    }

    if (rank == 0) std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    MPI_Bcast(&failures, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return failures == 0 ? 0 : 1;
}
