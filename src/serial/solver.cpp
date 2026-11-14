// Function bodies for serial/solver.hpp.

#include "serial/solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <optional>
#include <utility>

#include "common/affinity.hpp"
#include "common/corners.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/poisson.hpp"
#include "common/stability.hpp"
#include "common/stencil.hpp"
#include "serial/serial_flow.hpp"

SerialResult run_serial(const Config& cfg) {
    const double dx = cfg.lx / cfg.nx;
    const double dy = cfg.ly / cfg.ny;
    const double dt = time_step(cfg);
    const int spinup = spinup_steps(cfg, dt);
    const Block block = whole_domain(cfg);

    Field current(cfg.ny, cfg.nx);
    Field next(cfg.ny, cfg.nx);
    FaceVelocities vel(cfg.ny, cfg.nx);
    SerialFlow flow(cfg, block, dx, dy, dt);
    const StencilCoeffs k = make_stencil_coeffs(cfg.kappa, dx, dy, dt);

    double u_seen = 0.0;
    double v_seen = 0.0;
    auto check = [&](long s) {
        if (cfg.check_interval > 0 && s % cfg.check_interval == 0) {
            double u = 0.0;
            double v = 0.0;
            max_face_speeds(vel, u, v);
            u_seen = std::max(u_seen, u);
            v_seen = std::max(v_seen, v);
            check_time_step(cfg, dt, u, v);
        }
    };

    // Spin-up: the flow alone, from rest.
    const SteadyClock::time_point spinup_start = SteadyClock::now();
    for (long s = 0; s < spinup; ++s) {
        flow.step(s, vel);
        check(s + 1);
    }
    const double spinup_seconds = seconds_since(spinup_start);

    set_initial_condition(current, cfg, block);
    const double mass0 = local_mass(current, dx, dy);
    const double abs_mass0 = local_abs_mass(current, dx, dy);

    auto snapshot = [&](int n) {
        Snapshot shot;
        shot.step = n;
        shot.time = (spinup + n) * dt;
        shot.mass = local_mass(current, dx, dy);
        local_energy_enstrophy(flow.psi(), flow.zeta(), block, dx, dy, shot.energy, shot.enstrophy);
        shot.pigment = copy_interior(current);
        shot.psi = copy_corners(flow.psi());
        shot.zeta = copy_corners(flow.zeta());
        return shot;
    };

    std::optional<SnapshotWriter> snapshots;
    if (cfg.snapshot_interval > 0) {
        snapshots.emplace(cfg, dt, spinup * dt, serial_solver_name);
        snapshots->write(snapshot(0));
    }

    // The Poisson solves of the timed loop; the spin-up's are part of its time.
    double poisson_seconds = 0.0;
    long poisson_iterations = 0;
    int poisson_iterations_max = 0;
    auto flow_step = [&](long s) {
        const SerialPoissonWork work = flow.step(s, vel);
        poisson_seconds += work.seconds;
        poisson_iterations += work.iterations;
        poisson_iterations_max = std::max(poisson_iterations_max, work.iterations);
    };

    double snapshot_seconds = 0.0;
    const SteadyClock::time_point start = SteadyClock::now();
    for (int n = 0; n < cfg.steps; ++n) {
        // 1. Pigment: no-flux wall ghosts, then the update with the velocities of psi^n.
        fill_wall_ghosts(current, block);
        for (int i = 1; i <= cfg.ny; ++i) {
            update_row(current, next, vel, k, i, 1, cfg.nx);
        }
        std::swap(current, next);

        // 2.-4. The flow: zeta^{n+1}, psi^{n+1} and the new face velocities.
        flow_step(spinup + n);
        check(spinup + n + 1);

        if (snapshots && ((n + 1) % cfg.snapshot_interval == 0 || n + 1 == cfg.steps)) {
            const SteadyClock::time_point snapshot_start = SteadyClock::now();
            snapshots->write(snapshot(n + 1));
            snapshot_seconds += seconds_since(snapshot_start);
        }
    }
    const double loop_seconds = seconds_since(start) - snapshot_seconds;

    SerialResult result;
    RunSummary& s = result.summary;
    s.solver = serial_solver_name;
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
    s.poisson_solves = cfg.steps;
    s.snapshot_s = snapshot_seconds;
    s.mass_initial = mass0;
    s.mass_final = local_mass(current, dx, dy);
    s.mass_error = abs_mass0 > 0 ? std::abs(s.mass_final - mass0) / abs_mass0 : std::abs(s.mass_final - mass0);
    local_energy_enstrophy(flow.psi(), flow.zeta(), block, dx, dy, s.energy, s.enstrophy);
    s.u_max_measured = u_seen;
    s.v_max_measured = v_seen;
    s.peclet = peclet(cfg);
    s.peclet_grid = peclet_grid(cfg, std::max(u_seen, v_seen));
    if (cfg.initial_case == InitialCase::cosine_decay) {
        const CosineErrors e = cosine_decay_errors(current, cfg, block, dt, cfg.steps);
        s.cosine_discrete_error = e.discrete;
        s.cosine_continuous_error = e.continuous;
    }
    s.compiler = compiler_id();
    s.build_type = build_type();

    result.field = copy_interior(current);
    return result;
}
