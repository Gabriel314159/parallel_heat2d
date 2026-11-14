// Function bodies for parallel3d/solver.hpp.
//
// One time step on every rank (docs/design_3d.md, "Model A", "One time step"):
//   1. the pigment of every layer of the block, with the face velocities of
//      psi^n: wall ghosts on the six walls, the 3D halo exchange (blocking, or
//      overlapped with the update of the strict interior), the 3D update;
//   2.-4. the 2D flow of the rank's x-y block, exactly as in heat_parallel
//         (parallel/parallel_flow.hpp), unless --freeze-flow keeps the flow
//         of the end of the spin-up.

#include "parallel3d/solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#ifdef __linux__
#include <sys/prctl.h>
#endif

#include "common/affinity.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/stability.hpp"
#include "common3d/field3d.hpp"
#include "common3d/io3d.hpp"
#include "common3d/physics3d.hpp"
#include "common3d/stencil3d.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"
#include "parallel/parallel_flow.hpp"
#include "parallel/parallel_tools.hpp"
#include "parallel3d/gather3d.hpp"
#include "parallel3d/halo3d.hpp"

namespace {

const char* const program = "heat3d_parallel";

// Every layer of the block outside the strict interior, each cell exactly once.
void update_boundary_3d(const Field3D& c, Field3D& next, const FaceVelocities& vel, const StencilCoeffs3D& k,
                        const Block& b, const Interior& in) {
    if (in.empty()) {
        update_region_3d(c, next, vel, k, 1, b.ny, 1, b.nx);
        return;
    }
    update_region_3d(c, next, vel, k, 1, in.i_lo - 1, 1, b.nx);             // north strip
    update_region_3d(c, next, vel, k, in.i_hi + 1, b.ny, 1, b.nx);          // south strip
    update_region_3d(c, next, vel, k, in.i_lo, in.i_hi, 1, in.j_lo - 1);     // west strip
    update_region_3d(c, next, vel, k, in.i_lo, in.i_hi, in.j_hi + 1, b.nx);  // east strip
}

// This rank's x-y block, one layer deep or all nz layers deep.
Box block_box(const Block& b, int layers) {
    return Box{0, b.row_offset, b.col_offset, layers, b.ny, b.nx};
}

}  // namespace

Parallel3DResult run_parallel_3d(const Config3D& c3, MPI_Comm world) {
    const Config& cfg = c3.base;
    ProcessLayout layout = make_layout(cfg, world);
    const LayoutGuard guard(layout);
    const Block& b = layout.block;
    const int nz = c3.nz;

    // Before the first OpenMP region, so that all loops below use this count.
    const int threads = configure_threads(cfg.threads.value_or(0));
    check_same_threads(threads, layout.comm);
#if defined(__linux__) && defined(PR_SET_THP_DISABLE)
    // As in heat_parallel: 4 KiB pages keep the first-touch placement per thread.
    if (threads > 1) prctl(PR_SET_THP_DISABLE, 1, 0, 0, 0);
#endif
    const int nodes = count_nodes(layout);
    const int cores = min_over_ranks(allowed_cores(), layout.comm);
    const Binding binding = binding_over_nodes(layout, thread_cores());
    if (layout.rank == 0) warn_about_binding(program, cores, binding, threads);

    const double dx = cfg.lx / cfg.nx;
    const double dy = cfg.ly / cfg.ny;
    const double dz = dz_of(c3);
    const double dt = time_step_3d(c3);
    const int spinup = spinup_steps(cfg, dt);

    Field3D current(nz, b.ny, b.nx);
    Field3D next(nz, b.ny, b.nx);
    FaceVelocities vel(b.ny, b.nx);
    const StencilCoeffs3D k = make_stencil_coeffs_3d(cfg.kappa, dx, dy, dz, dt);
    HaloExchange halos(layout, b.ny, b.nx);  // the 2D flow's corner fields
    Halo3D pigment_halos(layout, nz, b.ny, b.nx);
    ParallelFlow flow(cfg, layout, dx, dy, dt, halos);
    const Interior interior = strict_interior(layout);

    // The section row, if this rank holds it (local index), else 0.
    const int section = section_row(c3);
    const int section_local =
        section > b.row_offset && section <= b.row_offset + b.ny ? section - b.row_offset : 0;

    // Every rank sees the same global maxima, so all ranks stop together.
    double u_seen = 0.0;
    double v_seen = 0.0;
    auto check_now = [&]() {
        double local[2] = {0.0, 0.0};
        max_face_speeds(vel, local[0], local[1]);
        double global[2] = {0.0, 0.0};
        mpi_check(MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, layout.comm), "speed check");
        u_seen = std::max(u_seen, global[0]);
        v_seen = std::max(v_seen, global[1]);
        check_time_step_3d(c3, dt, global[0], global[1]);
    };
    auto check = [&](long s) {
        if (cfg.check_interval > 0 && s % cfg.check_interval == 0) check_now();
    };

    // Spin-up: the flow alone, from rest.
    mpi_check(MPI_Barrier(layout.comm), "barrier before the spin-up");
    const double spinup_start = MPI_Wtime();
    for (long s = 0; s < spinup; ++s) {
        flow.step(s, vel);
        check(s + 1);
    }
    const double local_spinup = MPI_Wtime() - spinup_start;
    if (c3.freeze_flow) check_now();  // the flow the pigment will see for the whole run

    set_initial_condition_3d(current, c3, b);
    const double mass0 = sum_over_ranks(local_mass_3d(current, dx, dy, dz), layout.comm);
    const double abs_mass0 = sum_over_ranks(local_abs_mass_3d(current, dx, dy, dz), layout.comm);

    // Collective: every rank takes part in the sums and gathers; rank 0 writes.
    auto snapshot = [&](int n, std::vector<double>& section_values) {
        Snapshot shot;
        shot.step = n;
        shot.time = (spinup + n) * dt;
        shot.mass = sum_over_ranks(local_mass_3d(current, dx, dy, dz), layout.comm);
        double e = 0.0, z = 0.0;
        local_energy_enstrophy(flow.psi(), flow.zeta(), b, dx, dy, e, z);
        shot.energy = sum_over_ranks(e, layout.comm);
        shot.enstrophy = sum_over_ranks(z, layout.comm);
        shot.pigment = gather_boxes(layout.comm, block_box(b, 1), copy_layer(current, 1), 1, cfg.ny, cfg.nx);
        const Box section_box = section_local > 0 ? Box{0, 0, b.col_offset, nz, 1, b.nx} : Box{};
        section_values = gather_boxes(layout.comm, section_box,
                                      section_local > 0 ? copy_section(current, section_local) : std::vector<double>{},
                                      nz, 1, cfg.nx);
        const std::vector<double> psi = gather_boxes(layout.comm, block_box(b, 1), copy_interior(flow.psi()), 1, cfg.ny, cfg.nx);
        const std::vector<double> zeta =
            gather_boxes(layout.comm, block_box(b, 1), copy_interior(flow.zeta()), 1, cfg.ny, cfg.nx);
        if (layout.rank == 0) {
            shot.psi = corners_from_entries(psi, cfg.ny, cfg.nx);
            shot.zeta = corners_from_entries(zeta, cfg.ny, cfg.nx);
        }
        return shot;
    };

    std::optional<SnapshotWriter3D> snapshots;
    double snapshot_seconds = 0.0;
    if (cfg.snapshot_interval > 0) {
        std::vector<double> section_values;
        const Snapshot shot = snapshot(0, section_values);
        root_io(program, layout, [&] {
            snapshots.emplace(c3, dt, spinup * dt, program);
            snapshots->write(shot, section_values);
        });
    }

    // The Poisson solves of the timed loop; the spin-up's are part of its time.
    double poisson_seconds = 0.0;
    CommTimes poisson_comm;
    long poisson_iterations = 0;
    int poisson_iterations_max = 0;
    double pigment_seconds = 0.0;
    double pigment_halo_seconds = 0.0;

    mpi_check(MPI_Barrier(layout.comm), "barrier before the timed loop");
    const double start = MPI_Wtime();
    for (int n = 0; n < cfg.steps; ++n) {
        // 1. Pigment, every layer, with the velocities of psi^n.
        const double pigment_start = MPI_Wtime();
        if (cfg.exchange == Exchange::blocking) {
            fill_wall_ghosts_3d(current, b);
            const double halo_start = MPI_Wtime();
            pigment_halos.exchange(current);
            pigment_halo_seconds += MPI_Wtime() - halo_start;
            update_region_3d(current, next, vel, k, 1, b.ny, 1, b.nx);
        } else {
            // Messages travel while the wall ghosts are filled and the interior is updated.
            double halo_start = MPI_Wtime();
            pigment_halos.start(current);
            pigment_halo_seconds += MPI_Wtime() - halo_start;
            fill_wall_ghosts_3d(current, b);
            if (!interior.empty()) {
                update_region_3d(current, next, vel, k, interior.i_lo, interior.i_hi, interior.j_lo, interior.j_hi);
            }
            halo_start = MPI_Wtime();
            pigment_halos.finish(current);
            pigment_halo_seconds += MPI_Wtime() - halo_start;
            update_boundary_3d(current, next, vel, k, b, interior);
        }
        std::swap(current, next);
        pigment_seconds += MPI_Wtime() - pigment_start;

        // 2.-4. The 2D flow: zeta^{n+1}, psi^{n+1} and the new face velocities.
        if (!c3.freeze_flow) {
            const PoissonWork work = flow.step(spinup + n, vel);
            poisson_seconds += work.seconds;
            poisson_comm.halo += work.comm.halo;
            poisson_comm.sum += work.comm.sum;
            poisson_iterations += work.iterations;
            poisson_iterations_max = std::max(poisson_iterations_max, work.iterations);
            check(spinup + n + 1);
        }

        if (cfg.snapshot_interval > 0 && ((n + 1) % cfg.snapshot_interval == 0 || n + 1 == cfg.steps)) {
            const double snapshot_start = MPI_Wtime();
            std::vector<double> section_values;
            const Snapshot shot = snapshot(n + 1, section_values);
            root_io(program, layout, [&] { snapshots->write(shot, section_values); });
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
    const double pigment_max = max_over_ranks(pigment_seconds, layout.comm);
    const double pigment_halo_max = max_over_ranks(pigment_halo_seconds, layout.comm);
    const double mass_final = sum_over_ranks(local_mass_3d(current, dx, dy, dz), layout.comm);
    double e = 0.0, z = 0.0;
    local_energy_enstrophy(flow.psi(), flow.zeta(), b, dx, dy, e, z);
    const double energy = sum_over_ranks(e, layout.comm);
    const double enstrophy = sum_over_ranks(z, layout.comm);

    Parallel3DResult result;
    result.is_root = layout.rank == 0;
    RunSummary& s = result.summary;
    if (cfg.initial_case == InitialCase::cosine_decay) {
        const CosineErrors errors = cosine_decay_errors_3d(current, c3, b, dt, cfg.steps);
        s.cosine_discrete_error = max_over_ranks(errors.discrete, layout.comm);
        s.cosine_continuous_error = max_over_ranks(errors.continuous, layout.comm);
    }
    // The whole field only when it is written: benchmark runs gather nothing.
    if (!cfg.output.empty()) {
        result.field = gather_boxes(layout.comm, block_box(b, nz), copy_interior_3d(current), nz, cfg.ny, cfg.nx);
    }
    const std::string library = mpi_library();

    s.solver = program;
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
    s.poisson_iterations = poisson_iterations;
    s.poisson_iterations_max = poisson_iterations_max;
    s.poisson_solves = c3.freeze_flow ? 0 : cfg.steps;
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
    describe_3d(c3, Report3D{pigment_max, pigment_halo_max}, s);
    return result;
}
