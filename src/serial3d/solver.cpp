// Function bodies for serial3d/solver.hpp.
//
// One time step (docs/design_3d.md, "Model A", "One time step"):
//   1. the pigment, every layer, with the face velocities of psi^n;
//   2.-4. the 2D flow, exactly as in heat_serial (serial/serial_flow.hpp),
//         unless --freeze-flow keeps the flow of the end of the spin-up.

#include "serial3d/solver.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>

#include "common/affinity.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/stability.hpp"
#include "common3d/field3d.hpp"
#include "common3d/io3d.hpp"
#include "common3d/physics3d.hpp"
#include "common3d/stencil3d.hpp"
#include "serial/serial_flow.hpp"

Serial3DResult run_serial_3d(const Config3D& c3) {
    const Config& cfg = c3.base;
    const int nz = c3.nz;
    const double dx = cfg.lx / cfg.nx;
    const double dy = cfg.ly / cfg.ny;
    const double dz = dz_of(c3);
    const double dt = time_step_3d(c3);
    const int spinup = spinup_steps(cfg, dt);
    const Block block = whole_domain(cfg);
    const int section = section_row(c3);

    Field3D current(nz, cfg.ny, cfg.nx);
    Field3D next(nz, cfg.ny, cfg.nx);
    FaceVelocities vel(cfg.ny, cfg.nx);
    SerialFlow flow(cfg, block, dx, dy, dt);
    const StencilCoeffs3D k = make_stencil_coeffs_3d(cfg.kappa, dx, dy, dz, dt);

    double u_seen = 0.0;
    double v_seen = 0.0;
    auto check = [&](long s) {
        if (cfg.check_interval > 0 && s % cfg.check_interval == 0) {
            double u = 0.0;
            double v = 0.0;
            max_face_speeds(vel, u, v);
            u_seen = std::max(u_seen, u);
            v_seen = std::max(v_seen, v);
            check_time_step_3d(c3, dt, u, v);
        }
    };

    // Spin-up: the flow alone, from rest.
    const SteadyClock::time_point spinup_start = SteadyClock::now();
    for (long s = 0; s < spinup; ++s) {
        flow.step(s, vel);
        check(s + 1);
    }
    const double spinup_seconds = seconds_since(spinup_start);
    if (c3.freeze_flow) {
        // The flow the pigment will see for the whole run: check it once.
        double u = 0.0;
        double v = 0.0;
        max_face_speeds(vel, u, v);
        u_seen = std::max(u_seen, u);
        v_seen = std::max(v_seen, v);
        check_time_step_3d(c3, dt, u, v);
    }

    set_initial_condition_3d(current, c3, block);
    const double mass0 = local_mass_3d(current, dx, dy, dz);
    const double abs_mass0 = local_abs_mass_3d(current, dx, dy, dz);

    auto snapshot = [&](int n) {
        Snapshot shot;
        shot.step = n;
        shot.time = (spinup + n) * dt;
        shot.mass = local_mass_3d(current, dx, dy, dz);
        local_energy_enstrophy(flow.psi(), flow.zeta(), block, dx, dy, shot.energy, shot.enstrophy);
        shot.pigment = copy_layer(current, 1);
        shot.psi = copy_corners(flow.psi());
        shot.zeta = copy_corners(flow.zeta());
        return shot;
    };

    std::optional<SnapshotWriter3D> snapshots;
    if (cfg.snapshot_interval > 0) {
        snapshots.emplace(c3, dt, spinup * dt, serial3d_solver_name);
        snapshots->write(snapshot(0), copy_section(current, section));
    }

    // The Poisson solves of the timed loop; the spin-up's are part of its time.
    double poisson_seconds = 0.0;
    long poisson_iterations = 0;
    int poisson_iterations_max = 0;
    double pigment_seconds = 0.0;

    double snapshot_seconds = 0.0;
    const SteadyClock::time_point start = SteadyClock::now();
    for (int n = 0; n < cfg.steps; ++n) {
        // 1. Pigment: no-flux wall ghosts on all six walls, then every layer.
        const SteadyClock::time_point pigment_start = SteadyClock::now();
        fill_wall_ghosts_3d(current, block);
        update_region_3d(current, next, vel, k, 1, cfg.ny, 1, cfg.nx);
        std::swap(current, next);
        pigment_seconds += seconds_since(pigment_start);

        // 2.-4. The 2D flow: zeta^{n+1}, psi^{n+1} and the new face velocities.
        if (!c3.freeze_flow) {
            const SerialPoissonWork work = flow.step(spinup + n, vel);
            poisson_seconds += work.seconds;
            poisson_iterations += work.iterations;
            poisson_iterations_max = std::max(poisson_iterations_max, work.iterations);
            check(spinup + n + 1);
        }

        if (snapshots && ((n + 1) % cfg.snapshot_interval == 0 || n + 1 == cfg.steps)) {
            const SteadyClock::time_point snapshot_start = SteadyClock::now();
            snapshots->write(snapshot(n + 1), copy_section(current, section));
            snapshot_seconds += seconds_since(snapshot_start);
        }
    }
    const double loop_seconds = seconds_since(start) - snapshot_seconds;

    Serial3DResult result;
    RunSummary& s = result.summary;
    s.solver = serial3d_solver_name;
    s.dt = dt;
    s.spinup_steps = spinup;
    s.t_release = spinup * dt;
    s.t_end = (spinup + cfg.steps) * dt;
    s.elapsed_s = loop_seconds;
    s.spinup_s = spinup_seconds;
    s.cores_per_rank = allowed_cores();
    s.binding = binding_name(node_binding(thread_cores()));
    s.poisson_s = poisson_seconds;
    s.poisson_iterations = poisson_iterations;
    s.poisson_iterations_max = poisson_iterations_max;
    s.poisson_solves = c3.freeze_flow ? 0 : cfg.steps;
    s.snapshot_s = snapshot_seconds;
    s.mass_initial = mass0;
    s.mass_final = local_mass_3d(current, dx, dy, dz);
    s.mass_error = abs_mass0 > 0 ? std::abs(s.mass_final - mass0) / abs_mass0 : std::abs(s.mass_final - mass0);
    local_energy_enstrophy(flow.psi(), flow.zeta(), block, dx, dy, s.energy, s.enstrophy);
    s.u_max_measured = u_seen;
    s.v_max_measured = v_seen;
    s.peclet = peclet(cfg);
    s.peclet_grid = peclet_grid(cfg, std::max(u_seen, v_seen));
    if (cfg.initial_case == InitialCase::cosine_decay) {
        const CosineErrors e = cosine_decay_errors_3d(current, c3, block, dt, cfg.steps);
        s.cosine_discrete_error = e.discrete;
        s.cosine_continuous_error = e.continuous;
    }
    s.compiler = compiler_id();
    s.build_type = build_type();
    describe_3d(c3, Report3D{pigment_seconds, 0.0}, s);

    result.field = copy_interior_3d(current);
    return result;
}
