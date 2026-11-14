// The wind-driven flow of the whole grid, advanced one step at a time:
// shared by the 2D and 3D serial solvers (heat_serial and heat3d_serial). In
// the 3D pigment model the flow is the 2D flow of this file, the same in every
// layer.

#pragma once

#include <chrono>
#include <cstddef>
#include <utility>
#include <vector>

#include "common/config.hpp"
#include "common/corners.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"
#include "common/poisson.hpp"

using SteadyClock = std::chrono::steady_clock;

inline double seconds_since(SteadyClock::time_point start) {
    return std::chrono::duration<double>(SteadyClock::now() - start).count();
}

// Time and CG iterations of one Poisson solve.
struct SerialPoissonWork {
    double seconds = 0.0;
    int iterations = 0;
};

// The flow state of the whole grid and one time step of it
// (docs/design.md, "One time step", steps 2 to 4).
class SerialFlow {
public:
    SerialFlow(const Config& cfg, const Block& block, double dx, double dy, double dt)
        : cfg_(cfg), block_(block), dx_(dx), dy_(dy), dt_(dt),
          psi_(block.ny, block.nx), zeta_(block.ny, block.nx),
          r0_(block.ny, block.nx), r1_(block.ny, block.nx), r2_(block.ny, block.nx),
          coeffs_(make_flow_coeffs(cfg.viscosity, dx, dy)),
          profile_(wind_profile(cfg, block)), unknown_(unknown_corners(block)),
          poisson_(block, dx, dy, cfg.poisson_tol, cfg.poisson_max_iter) {}

    // Advance zeta and psi from flow step s (time s dt) to s + 1, and refresh
    // the face velocities. Returns the time and iterations of the Poisson solve.
    SerialPoissonWork step(long s, FaceVelocities& vel) {
        const double strength = wind_strength(cfg_, s * dt_);
        for (int i = unknown_.i_first; i <= unknown_.i_last; ++i) {
            vorticity_rhs_row(psi_, zeta_, r0_, coeffs_, strength * profile_[static_cast<std::size_t>(i)], nullptr,
                              i, unknown_.j_first, unknown_.j_last);
        }
        const AbWeights w = ab_weights(s);
        for (int i = unknown_.i_first; i <= unknown_.i_last; ++i) {
            ab_update_row(zeta_, r0_, r1_, r2_, w, dt_, i, unknown_.j_first, unknown_.j_last);
        }
        // psi still holds psi^n: the warm start of the solve.
        const SteadyClock::time_point start = SteadyClock::now();
        const CgStats stats = poisson_.solve(zeta_, psi_, comm_);
        const SerialPoissonWork work{seconds_since(start), stats.iterations};
        faces_from_psi(psi_, vel, dx_, dy_);
        // The newest right-hand side becomes R^{n-1} of the next step.
        std::swap(r2_, r1_);
        std::swap(r1_, r0_);
        return work;
    }

    const Field& psi() const { return psi_; }
    const Field& zeta() const { return zeta_; }

private:
    const Config& cfg_;
    Block block_;
    double dx_;
    double dy_;
    double dt_;
    Field psi_;   // corner fields; the wall corners and ghosts stay 0
    Field zeta_;
    Field r0_;    // R^n, R^{n-1}, R^{n-2}
    Field r1_;
    Field r2_;
    FlowCoeffs coeffs_;
    std::vector<double> profile_;
    CornerRange unknown_;
    PoissonCG poisson_;
    SerialComm comm_;
};
