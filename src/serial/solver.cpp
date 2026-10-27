// Function bodies for serial/solver.hpp.

#include "serial/solver.hpp"

#include <chrono>
#include <cmath>
#include <optional>
#include <utility>

#include "common/grid.hpp"
#include "common/physics.hpp"
#include "common/stencil.hpp"
#include "common/velocity.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

}  // namespace

SerialResult run_serial(const Config& cfg) {
    const double dx = cfg.lx / cfg.nx;
    const double dy = cfg.ly / cfg.ny;
    const double dt = time_step(cfg);
    const Block block = whole_domain(cfg);

    Field current(cfg.ny, cfg.nx);
    Field next(cfg.ny, cfg.nx);
    FaceVelocities vel(cfg.ny, cfg.nx);
    FlowField flow(flow_params(cfg), block, dx, dy);
    const StencilCoeffs k = make_stencil_coeffs(cfg.kappa, dx, dy, dt);

    set_initial_condition(current, cfg, block);
    const double mass0 = local_mass(current, dx, dy);
    const double abs_mass0 = local_abs_mass(current, dx, dy);

    // A steady flow needs its face velocities only once.
    if (flow.steady()) flow.fill_all(vel);

    std::optional<SnapshotWriter> snapshots;
    if (cfg.snapshot_interval > 0) {
        snapshots.emplace(cfg, dt, "heat_serial");
        snapshots->write(0, copy_interior(current), mass0);
    }

    double snapshot_seconds = 0.0;
    const Clock::time_point start = Clock::now();
    for (int n = 0; n < cfg.steps; ++n) {
        // 1. No-flux walls: refill the ghost cells from the current interior.
        fill_wall_ghosts(current, block);

        // 2. Time-dependent flow: face velocities at the old time t_n.
        if (!flow.steady()) {
            flow.set_time(n * dt);
            flow.fill_all(vel);
        }

        // 3. Update every interior cell.
        for (int i = 1; i <= cfg.ny; ++i) {
            update_row(current, next, vel, k, i, 1, cfg.nx);
        }

        // 4. The new field becomes the current one (no copying, see grid.hpp).
        std::swap(current, next);

        if (snapshots && ((n + 1) % cfg.snapshot_interval == 0 || n + 1 == cfg.steps)) {
            const Clock::time_point snapshot_start = Clock::now();
            snapshots->write(n + 1, copy_interior(current), local_mass(current, dx, dy));
            snapshot_seconds += seconds_since(snapshot_start);
        }
    }
    const double loop_seconds = seconds_since(start) - snapshot_seconds;

    SerialResult result;
    RunSummary& s = result.summary;
    s.solver = "heat_serial";
    s.dt = dt;
    s.t_end = cfg.steps * dt;
    s.elapsed_s = loop_seconds;
    s.snapshot_s = snapshot_seconds;
    s.mass_initial = mass0;
    s.mass_final = local_mass(current, dx, dy);
    s.mass_error = abs_mass0 > 0 ? std::abs(s.mass_final - mass0) / abs_mass0 : std::abs(s.mass_final - mass0);
    s.peclet = peclet(cfg);
    s.peclet_grid = peclet_grid(cfg);
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
