// Function bodies for parallel/solver.hpp.

#include "parallel/solver.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

#ifdef HEAT2D_USE_OPENMP
#include <omp.h>
#endif

#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/stencil.hpp"
#include "common/velocity.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"

namespace {

// Update rows i0..i1, columns j0..j1 (inclusive), rows shared out among the
// OpenMP threads. Each destination cell has exactly one writer.
void update_region(const Field& c, Field& next, const FaceVelocities& vel, const StencilCoeffs& k,
                   int i0, int i1, int j0, int j1) {
    if (i0 > i1 || j0 > j1) return;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = i0; i <= i1; ++i) {
        update_row(c, next, vel, k, i, j0, j1);
    }
}

// Face velocities of the block at the time of the last flow.set_time().
void fill_velocities(const FlowField& flow, FaceVelocities& vel, int ny) {
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 1; i <= ny + 1; ++i) {
        if (i <= ny) flow.fill_east_row(vel, i);
        flow.fill_north_row(vel, i);
    }
}

// The strict interior for overlap (src/README.md, "Strict interior under
// overlap"): cells at least `halo` away from every edge shared with another
// rank. Cells next to a physical wall can be interior, because their ghosts
// are filled locally before the update.
struct Interior {
    int i_lo;
    int i_hi;
    int j_lo;
    int j_hi;
    bool empty() const { return i_lo > i_hi || j_lo > j_hi; }
};

Interior strict_interior(const ProcessLayout& layout) {
    const Block& b = layout.block;
    return Interior{layout.north != MPI_PROC_NULL ? 1 + halo : 1, layout.south != MPI_PROC_NULL ? b.ny - halo : b.ny,
                    layout.west != MPI_PROC_NULL ? 1 + halo : 1, layout.east != MPI_PROC_NULL ? b.nx - halo : b.nx};
}

// Everything outside the interior, each cell exactly once.
void update_boundary(const Field& c, Field& next, const FaceVelocities& vel, const StencilCoeffs& k,
                     const Block& b, const Interior& in) {
    if (in.empty()) {
        update_region(c, next, vel, k, 1, b.ny, 1, b.nx);
        return;
    }
    update_region(c, next, vel, k, 1, in.i_lo - 1, 1, b.nx);         // top strip
    update_region(c, next, vel, k, in.i_hi + 1, b.ny, 1, b.nx);      // bottom strip
    update_region(c, next, vel, k, in.i_lo, in.i_hi, 1, in.j_lo - 1);  // left strip
    update_region(c, next, vel, k, in.i_lo, in.i_hi, in.j_hi + 1, b.nx);  // right strip
}

// Collects every rank's interior into one global field on rank 0.
class FieldGatherer {
public:
    FieldGatherer(const ProcessLayout& layout, int ny, int nx)
        : comm_(layout.comm), rank_(layout.rank), size_(layout.size), ny_(ny), nx_(nx) {
        const Block& b = layout.block;
        const int mine[4] = {b.row_offset, b.col_offset, b.ny, b.nx};
        if (rank_ == 0) blocks_.resize(static_cast<std::size_t>(4) * size_);
        MPI_Gather(mine, 4, MPI_INT, rank_ == 0 ? blocks_.data() : nullptr, 4, MPI_INT, 0, comm_);
        if (rank_ == 0) {
            counts_.resize(size_);
            displacements_.resize(size_);
            int offset = 0;
            for (int r = 0; r < size_; ++r) {
                counts_[r] = blocks_[4 * r + 2] * blocks_[4 * r + 3];
                displacements_[r] = offset;
                offset += counts_[r];
            }
        }
    }

    // Collective. The global ny x nx field on rank 0, an empty vector elsewhere.
    std::vector<double> gather(const Field& c) const {
        const std::vector<double> local = copy_interior(c);
        std::vector<double> received;
        if (rank_ == 0) received.resize(static_cast<std::size_t>(ny_) * nx_);
        MPI_Gatherv(local.data(), static_cast<int>(local.size()), MPI_DOUBLE, received.data(),
                    rank_ == 0 ? counts_.data() : nullptr, rank_ == 0 ? displacements_.data() : nullptr, MPI_DOUBLE,
                    0, comm_);
        if (rank_ != 0) return {};

        std::vector<double> global(static_cast<std::size_t>(ny_) * nx_);
        for (int r = 0; r < size_; ++r) {
            const int row0 = blocks_[4 * r];
            const int col0 = blocks_[4 * r + 1];
            const int rows = blocks_[4 * r + 2];
            const int cols = blocks_[4 * r + 3];
            const double* block = received.data() + displacements_[r];
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    global[static_cast<std::size_t>(row0 + i) * nx_ + (col0 + j)] =
                        block[static_cast<std::size_t>(i) * cols + j];
                }
            }
        }
        return global;
    }

private:
    MPI_Comm comm_;
    int rank_;
    int size_;
    int ny_;
    int nx_;
    std::vector<int> blocks_;  // per rank: row offset, column offset, rows, columns (rank 0 only)
    std::vector<int> counts_;
    std::vector<int> displacements_;
};

double sum_over_ranks(double value, MPI_Comm comm) {
    double total = 0.0;
    MPI_Allreduce(&value, &total, 1, MPI_DOUBLE, MPI_SUM, comm);
    return total;
}

double max_over_ranks(double value, MPI_Comm comm) {
    double result = 0.0;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
    return result;
}

// File output happens on rank 0 only. If it fails, the other ranks would wait
// forever in the next collective call, so stop the whole job instead.
template <typename Action>
void root_io(const ProcessLayout& layout, Action action) {
    if (layout.rank != 0) return;
    try {
        action();
    } catch (const std::exception& e) {
        std::cerr << "heat_parallel: error: " << e.what() << "\n";
        MPI_Abort(layout.comm, 3);
    }
}

int threads_in_use() {
    int threads = 1;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel
    {
#pragma omp single
        threads = omp_get_num_threads();
    }
#endif
    return threads;
}

std::string mpi_library() {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int length = 0;
    MPI_Get_library_version(version, &length);
    std::string text(version, static_cast<std::size_t>(length));
    const std::size_t end = text.find_first_of("\n,");
    return end == std::string::npos ? text : text.substr(0, end);
}

}  // namespace

ParallelResult run_parallel(const Config& cfg, MPI_Comm world) {
    ProcessLayout layout = make_layout(cfg, world);
    const LayoutGuard guard(layout);
    const Block& b = layout.block;

    const double dx = cfg.lx / cfg.nx;
    const double dy = cfg.ly / cfg.ny;
    const double dt = time_step(cfg);

    Field current(b.ny, b.nx);
    Field next(b.ny, b.nx);
    FaceVelocities vel(b.ny, b.nx);
    FlowField flow(flow_params(cfg), b, dx, dy);
    const StencilCoeffs k = make_stencil_coeffs(cfg.kappa, dx, dy, dt);
    HaloExchange halos(layout, b.ny, b.nx);
    const FieldGatherer gatherer(layout, cfg.ny, cfg.nx);
    const Interior interior = strict_interior(layout);

    set_initial_condition(current, cfg, b);
    const double mass0 = sum_over_ranks(local_mass(current, dx, dy), layout.comm);
    const double abs_mass0 = sum_over_ranks(local_abs_mass(current, dx, dy), layout.comm);

    if (flow.steady()) fill_velocities(flow, vel, b.ny);

    std::optional<SnapshotWriter> snapshots;
    double snapshot_seconds = 0.0;
    if (cfg.snapshot_interval > 0) {
        const std::vector<double> field = gatherer.gather(current);
        root_io(layout, [&] {
            snapshots.emplace(cfg, dt, "heat_parallel");
            snapshots->write(0, field, mass0);
        });
    }

    MPI_Barrier(layout.comm);
    const double start = MPI_Wtime();
    for (int n = 0; n < cfg.steps; ++n) {
        if (cfg.exchange == Exchange::blocking) {
            fill_wall_ghosts(current, b);
            halos.exchange(current);
            if (!flow.steady()) {
                flow.set_time(n * dt);
                fill_velocities(flow, vel, b.ny);
            }
            update_region(current, next, vel, k, 1, b.ny, 1, b.nx);
        } else {
            // Messages travel while the local work and the interior update run.
            halos.start(current);
            fill_wall_ghosts(current, b);
            if (!flow.steady()) {
                flow.set_time(n * dt);
                fill_velocities(flow, vel, b.ny);
            }
            if (!interior.empty()) {
                update_region(current, next, vel, k, interior.i_lo, interior.i_hi, interior.j_lo, interior.j_hi);
            }
            halos.finish(current);
            update_boundary(current, next, vel, k, b, interior);
        }
        std::swap(current, next);

        if (cfg.snapshot_interval > 0 && ((n + 1) % cfg.snapshot_interval == 0 || n + 1 == cfg.steps)) {
            const double snapshot_start = MPI_Wtime();
            const double mass = sum_over_ranks(local_mass(current, dx, dy), layout.comm);
            const std::vector<double> field = gatherer.gather(current);
            root_io(layout, [&] { snapshots->write(n + 1, field, mass); });
            snapshot_seconds += MPI_Wtime() - snapshot_start;
        }
    }
    const double local_loop = MPI_Wtime() - start - snapshot_seconds;

    // Everything below is outside the timed loop.
    double loop_seconds = 0.0;
    MPI_Reduce(&local_loop, &loop_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, layout.comm);
    const double mass_final = sum_over_ranks(local_mass(current, dx, dy), layout.comm);

    ParallelResult result;
    result.is_root = layout.rank == 0;
    RunSummary& s = result.summary;
    if (cfg.initial_case == InitialCase::cosine_decay) {
        const CosineErrors e = cosine_decay_errors(current, cfg, b, dt, cfg.steps);
        s.cosine_discrete_error = max_over_ranks(e.discrete, layout.comm);
        s.cosine_continuous_error = max_over_ranks(e.continuous, layout.comm);
    }
    result.field = gatherer.gather(current);
    const int threads = threads_in_use();
    const std::string library = mpi_library();

    s.solver = "heat_parallel";
    s.mpi_ranks = layout.size;
    s.omp_threads = threads;
    s.py = layout.py;
    s.px = layout.px;
    s.decomposition = decomp_name(cfg.decomposition);
    s.exchange = exchange_name(cfg.exchange);
    s.dt = dt;
    s.t_end = cfg.steps * dt;
    s.elapsed_s = loop_seconds;
    s.snapshot_s = snapshot_seconds;
    s.mass_initial = mass0;
    s.mass_final = mass_final;
    s.mass_error = abs_mass0 > 0 ? std::abs(mass_final - mass0) / abs_mass0 : std::abs(mass_final - mass0);
    s.peclet = peclet(cfg);
    s.peclet_grid = peclet_grid(cfg);
    s.compiler = compiler_id();
    s.mpi_library = library;
    s.build_type = build_type();
    return result;
}
