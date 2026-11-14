// Function bodies for parallel/solver.hpp.

#include "parallel/solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef HEAT2D_USE_OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sys/prctl.h>
#endif

#include "common/affinity.hpp"
#include "common/corners.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/stability.hpp"
#include "common/stencil.hpp"
#include "common/poisson.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"
#include "parallel/parallel_flow.hpp"
#include "parallel/parallel_tools.hpp"
#include "parallel/poisson_comm.hpp"

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
        mpi_check(MPI_Gather(mine, 4, MPI_INT, rank_ == 0 ? blocks_.data() : nullptr, 4, MPI_INT, 0, comm_),
                  "gather block sizes");
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
        mpi_check(MPI_Gatherv(local.data(), static_cast<int>(local.size()), MPI_DOUBLE, received.data(),
                              rank_ == 0 ? counts_.data() : nullptr, rank_ == 0 ? displacements_.data() : nullptr,
                              MPI_DOUBLE, 0, comm_),
                  "gather field");
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

}  // namespace

ParallelResult run_parallel(const Config& cfg, MPI_Comm world) {
    ProcessLayout layout = make_layout(cfg, world);
    const LayoutGuard guard(layout);
    const Block& b = layout.block;

#ifdef HEAT2D_USE_OFFLOAD
    // One GPU per rank: the ranks of a node take its GPUs in turn
    // (docs/design.md, "The Poisson solve on a GPU").
    if (cfg.exchange == Exchange::overlap) {
        throw std::invalid_argument("the GPU build supports only --exchange blocking");
    }
    const int devices = omp_get_num_devices();
    const int node_rank = rank_on_node(layout);
    if (devices > 0) omp_set_default_device(node_rank % devices);
    if (layout.rank == 0 && devices == 0) {
        std::cerr << "heat_parallel: warning: GPU build, but no GPU found; the Poisson solve runs on the host\n";
    }
#endif
    // Before the first OpenMP region, so that all loops below use this count.
    const int threads = configure_threads(cfg.threads.value_or(0));
    check_same_threads(threads, layout.comm);
#if defined(__linux__) && defined(PR_SET_THP_DISABLE)
    // With 2 MiB pages, whichever thread touches a page first would decide the
    // socket of a whole block of rows; 4 KiB pages keep the placement per thread.
    if (threads > 1) prctl(PR_SET_THP_DISABLE, 1, 0, 0, 0);
#endif
    const int nodes = count_nodes(layout);

    // Threads that share a core, or that the operating system may move from
    // core to core, make the timings meaningless per core: warn, but run, so
    // that the run report shows what happened.
    const int cores = min_over_ranks(allowed_cores(), layout.comm);
    const Binding binding = binding_over_nodes(layout, thread_cores());
    if (layout.rank == 0) warn_about_binding("heat_parallel", cores, binding, threads);

    const double dx = cfg.lx / cfg.nx;
    const double dy = cfg.ly / cfg.ny;
    const double dt = time_step(cfg);
    const int spinup = spinup_steps(cfg, dt);

    Field current(b.ny, b.nx);
    Field next(b.ny, b.nx);
    FaceVelocities vel(b.ny, b.nx);
    const StencilCoeffs k = make_stencil_coeffs(cfg.kappa, dx, dy, dt);
    HaloExchange halos(layout, b.ny, b.nx);
    ParallelFlow flow(cfg, layout, dx, dy, dt, halos);
    const FieldGatherer gatherer(layout, cfg.ny, cfg.nx);
    const Interior interior = strict_interior(layout);

    // Every rank sees the same global maxima, so all ranks stop together.
    double u_seen = 0.0;
    double v_seen = 0.0;
    auto check = [&](long s) {
        if (cfg.check_interval > 0 && s % cfg.check_interval == 0) {
            double local[2] = {0.0, 0.0};
            max_face_speeds(vel, local[0], local[1]);
            double global[2] = {0.0, 0.0};
            mpi_check(MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, layout.comm), "speed check");
            u_seen = std::max(u_seen, global[0]);
            v_seen = std::max(v_seen, global[1]);
            check_time_step(cfg, dt, global[0], global[1]);
        }
    };

    // Spin-up: the flow alone, from rest.
    mpi_check(MPI_Barrier(layout.comm), "barrier before the spin-up");
    const double spinup_start = MPI_Wtime();
    for (long s = 0; s < spinup; ++s) {
        flow.step(s, vel);
        check(s + 1);
    }
    const double local_spinup = MPI_Wtime() - spinup_start;

    set_initial_condition(current, cfg, b);
    const double mass0 = sum_over_ranks(local_mass(current, dx, dy), layout.comm);
    const double abs_mass0 = sum_over_ranks(local_abs_mass(current, dx, dy), layout.comm);

    // Collective: every rank takes part in the sums and gathers; rank 0 writes.
    auto snapshot = [&](int n) {
        Snapshot shot;
        shot.step = n;
        shot.time = (spinup + n) * dt;
        shot.mass = sum_over_ranks(local_mass(current, dx, dy), layout.comm);
        double e = 0.0, z = 0.0;
        local_energy_enstrophy(flow.psi(), flow.zeta(), b, dx, dy, e, z);
        shot.energy = sum_over_ranks(e, layout.comm);
        shot.enstrophy = sum_over_ranks(z, layout.comm);
        shot.pigment = gatherer.gather(current);
        const std::vector<double> psi = gatherer.gather(flow.psi());
        const std::vector<double> zeta = gatherer.gather(flow.zeta());
        if (layout.rank == 0) {
            shot.psi = corners_from_entries(psi, cfg.ny, cfg.nx);
            shot.zeta = corners_from_entries(zeta, cfg.ny, cfg.nx);
        }
        return shot;
    };

    std::optional<SnapshotWriter> snapshots;
    double snapshot_seconds = 0.0;
    if (cfg.snapshot_interval > 0) {
        const Snapshot shot = snapshot(0);
        root_io("heat_parallel", layout, [&] {
            snapshots.emplace(cfg, dt, spinup * dt, "heat_parallel");
            snapshots->write(shot);
        });
    }

    // The Poisson solves of the timed loop; the spin-up's are part of its time.
    double poisson_seconds = 0.0;
    double poisson_copy = 0.0;
    bool poisson_on_device = false;
    CommTimes poisson_comm;
    long poisson_iterations = 0;
    int poisson_iterations_max = 0;
    auto flow_step = [&](long s) {
        const PoissonWork work = flow.step(s, vel);
        poisson_seconds += work.seconds;
        poisson_comm.halo += work.comm.halo;
        poisson_comm.sum += work.comm.sum;
        poisson_copy += work.copy_seconds;
        poisson_on_device = work.on_device;
        poisson_iterations += work.iterations;
        poisson_iterations_max = std::max(poisson_iterations_max, work.iterations);
    };

    mpi_check(MPI_Barrier(layout.comm), "barrier before the timed loop");
    const double start = MPI_Wtime();
    for (int n = 0; n < cfg.steps; ++n) {
        // 1. Pigment with the velocities of psi^n.
        if (cfg.exchange == Exchange::blocking) {
            fill_wall_ghosts(current, b);
            halos.exchange(current);
            update_region(current, next, vel, k, 1, b.ny, 1, b.nx);
        } else {
            // Messages travel while the wall ghosts are filled and the interior is updated.
            halos.start(current);
            fill_wall_ghosts(current, b);
            if (!interior.empty()) {
                update_region(current, next, vel, k, interior.i_lo, interior.i_hi, interior.j_lo, interior.j_hi);
            }
            halos.finish(current);
            update_boundary(current, next, vel, k, b, interior);
        }
        std::swap(current, next);

        // 2.-4. The flow: zeta^{n+1}, psi^{n+1} and the new face velocities.
        flow_step(spinup + n);
        check(spinup + n + 1);

        if (cfg.snapshot_interval > 0 && ((n + 1) % cfg.snapshot_interval == 0 || n + 1 == cfg.steps)) {
            const double snapshot_start = MPI_Wtime();
            const Snapshot shot = snapshot(n + 1);
            root_io("heat_parallel", layout, [&] { snapshots->write(shot); });
            snapshot_seconds += MPI_Wtime() - snapshot_start;
        }
    }
    const double local_loop = MPI_Wtime() - start - snapshot_seconds;

    // Everything below is outside the timed loops.
    const double loop_seconds = max_over_ranks(local_loop, layout.comm);
    const double spinup_seconds = max_over_ranks(local_spinup, layout.comm);
    const double poisson_max = max_over_ranks(poisson_seconds, layout.comm);
    const double poisson_halo_max = max_over_ranks(poisson_comm.halo, layout.comm);
    const double poisson_sum_max = max_over_ranks(poisson_comm.sum, layout.comm);
    const double mass_final = sum_over_ranks(local_mass(current, dx, dy), layout.comm);
    double e = 0.0, z = 0.0;
    local_energy_enstrophy(flow.psi(), flow.zeta(), b, dx, dy, e, z);
    const double energy = sum_over_ranks(e, layout.comm);
    const double enstrophy = sum_over_ranks(z, layout.comm);

    ParallelResult result;
    result.is_root = layout.rank == 0;
    RunSummary& s = result.summary;
    if (cfg.initial_case == InitialCase::cosine_decay) {
        const CosineErrors errors = cosine_decay_errors(current, cfg, b, dt, cfg.steps);
        s.cosine_discrete_error = max_over_ranks(errors.discrete, layout.comm);
        s.cosine_continuous_error = max_over_ranks(errors.continuous, layout.comm);
    }
    result.field = gatherer.gather(current);
    const std::string library = mpi_library();

    s.solver = "heat_parallel";
    s.mpi_ranks = layout.size;
    s.omp_threads = threads;
    s.nodes = nodes;
    s.cores_per_rank = cores;
    s.binding = binding_name(binding);
    s.py = layout.py;
    s.px = layout.px;
    s.decomposition = decomp_name(cfg.decomposition);
    s.exchange = exchange_name(cfg.exchange);
    s.dt = dt;
    s.spinup_steps = spinup;
    s.t_release = spinup * dt;
    s.t_end = (spinup + cfg.steps) * dt;
    s.elapsed_s = loop_seconds;
    s.spinup_s = spinup_seconds;
    s.poisson_s = poisson_max;
    s.poisson_halo_s = poisson_halo_max;
    s.poisson_sum_s = poisson_sum_max;
    s.poisson_copy_s = max_over_ranks(poisson_copy, layout.comm);
#ifdef HEAT2D_USE_OFFLOAD
    // Every rank must have run its solve on a GPU for the run to count as one.
    const int on_device = static_cast<int>(min_over_ranks(poisson_on_device ? 1 : 0, layout.comm));
    s.poisson_device = cfg.steps == 0 ? "gpu (no solve)" : (on_device ? "gpu" : "host (no GPU found)");
#endif
    s.poisson_iterations = poisson_iterations;
    s.poisson_iterations_max = poisson_iterations_max;
    s.poisson_solves = cfg.steps;
    s.snapshot_s = snapshot_seconds;
    s.mass_initial = mass0;
    s.mass_final = mass_final;
    s.mass_error = abs_mass0 > 0 ? std::abs(mass_final - mass0) / abs_mass0 : std::abs(mass_final - mass0);
    s.energy = energy;
    s.enstrophy = enstrophy;
    s.u_max_measured = u_seen;
    s.v_max_measured = v_seen;
    s.peclet = peclet(cfg);
    s.peclet_grid = peclet_grid(cfg, std::max(u_seen, v_seen));
    s.compiler = compiler_id();
    s.mpi_library = library;
    s.build_type = build_type();
    return result;
}
