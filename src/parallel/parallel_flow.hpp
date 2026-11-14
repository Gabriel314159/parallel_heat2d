// The wind-driven flow of one rank's block, advanced one step at a time:
// shared by the 2D and 3D parallel solvers (heat_parallel and
// heat3d_parallel). In the 3D pigment model the flow is the 2D flow of this
// file, the same in every layer.
//
// Header-only: each solver compiles its loops with its own OpenMP setting.

#pragma once

#include <mpi.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "common/config.hpp"
#include "common/corners.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"
#include "common/poisson.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"
#include "parallel/poisson_comm.hpp"

// Time and CG iterations of one Poisson solve, and the time inside its
// communication (see CommTimes).
struct PoissonWork {
    double seconds = 0.0;
    int iterations = 0;
    CommTimes comm;
    double copy_seconds = 0.0;  // GPU build: copies between host and GPU
    bool on_device = false;     // GPU build: the solve ran on a GPU
};

// The flow state of one rank's block and one time step of it
// (docs/design.md, "One time step", steps 2 to 4), rows shared out among the
// OpenMP threads.
class ParallelFlow {
public:
    ParallelFlow(const Config& cfg, const ProcessLayout& layout, double dx, double dy, double dt, HaloExchange& halos)
        : cfg_(cfg), block_(layout.block), dx_(dx), dy_(dy), dt_(dt), halos_(halos),
          psi_(block_.ny, block_.nx), zeta_(block_.ny, block_.nx),
          r0_(block_.ny, block_.nx), r1_(block_.ny, block_.nx), r2_(block_.ny, block_.nx),
          coeffs_(make_flow_coeffs(cfg.viscosity, dx, dy)),
          profile_(wind_profile(cfg, block_)), unknown_(unknown_corners(block_)),
          poisson_(block_, dx, dy, cfg.poisson_tol, cfg.poisson_max_iter),
          comm_(layout.comm, halos, cfg.exchange == Exchange::overlap) {}

    // Advance zeta and psi from flow step s (time s dt) to s + 1, and refresh
    // the face velocities. Returns the time and iterations of the Poisson solve.
    PoissonWork step(long s, FaceVelocities& vel) {
        // zeta^n with its corner ghosts; psi^n already has them from the previous step.
        zero_wall_ghost_corners(zeta_, block_);
        halos_.exchange_with_corners(zeta_);

        const double strength = wind_strength(cfg_, s * dt_);
        const AbWeights w = ab_weights(s);
        const CornerRange u = unknown_;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int i = u.i_first; i <= u.i_last; ++i) {
            vorticity_rhs_row(psi_, zeta_, r0_, coeffs_, strength * profile_[static_cast<std::size_t>(i)], nullptr, i,
                              u.j_first, u.j_last);
        }
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int i = u.i_first; i <= u.i_last; ++i) {
            ab_update_row(zeta_, r0_, r1_, r2_, w, dt_, i, u.j_first, u.j_last);
        }

        // psi still holds psi^n: the warm start of the solve. All ranks take
        // the same number of iterations, because they see the same global sums.
        comm_.reset_times();
        const double start = MPI_Wtime();
        const CgStats stats = poisson_.solve(zeta_, psi_, comm_);
        const PoissonWork work{MPI_Wtime() - start, stats.iterations, comm_.times(), stats.copy_seconds, stats.on_device};

        zero_wall_ghost_corners(psi_, block_);
        halos_.exchange_with_corners(psi_);
        const int ny = block_.ny;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int i = 1; i <= ny + 1; ++i) {
            faces_from_psi_row(psi_, vel, dx_, dy_, i);
        }
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
    HaloExchange& halos_;
    Field psi_;
    Field zeta_;
    Field r0_;
    Field r1_;
    Field r2_;
    FlowCoeffs coeffs_;
    std::vector<double> profile_;
    CornerRange unknown_;
    PoissonCG poisson_;
    MpiComm comm_;
};
