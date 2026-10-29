// The Poisson equation for psi on the corners, solved with the conjugate
// gradient (CG) method. See docs/design.md, "Poisson equation for ψ".
//
// At every unknown corner,
//   A psi = -(psi_E - 2 psi + psi_W)/dx^2 - (psi_N - 2 psi + psi_S)/dy^2 = zeta,
// with psi = 0 at the wall corners. A is symmetric and positive definite and is
// never stored: applying it is the stencil above.
//
// The solver contains no MPI. What it needs from outside, a global sum and a
// halo exchange (blocking, or started and finished around the work that reads
// no halo), is passed in as a PoissonComm, so the serial and the parallel
// solver share this one implementation. It is header-only, so the parallel
// solver compiles its loops with OpenMP (HEAT2D_USE_OPENMP) and the serial
// solver without.

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "common/corners.hpp"
#include "common/grid.hpp"

// What the solver needs from its caller.
class PoissonComm {
public:
    virtual ~PoissonComm() = default;
    // Replace values[0..count) by their sums over all ranks.
    virtual void sum(double* values, int count) const = 0;
    // Fill the process halos of a corner field (only the four nearest
    // neighbours are read by A).
    virtual void exchange(Field& f) const = 0;
    // Overlap (--exchange overlap): start() begins filling the process halos
    // and returns at once, finish() waits until they are filled. In between,
    // the solver applies A where no process halo is read, and neither changes
    // the field nor reads its process halos. By default there is no overlap.
    virtual bool overlap() const { return false; }
    virtual void start(Field& f) const { exchange(f); }
    virtual void finish(Field&) const {}
};

// One process: local sums are global sums, and there are no process halos.
class SerialComm : public PoissonComm {
public:
    void sum(double*, int) const override {}
    void exchange(Field&) const override {}
};

struct CgStats {
    int iterations = 0;
    double relative_residual = 0.0;  // ||r|| / ||zeta|| at the end
};

// The unknown corners of a block, split for overlap. A reads one corner in
// each direction, so the `inner` corners, one corner away from every block
// edge shared with another rank (every side that is not a wall), read no
// process halo. The `edges` are the rest, each corner exactly once; some of
// them may be empty. A block too small to have inner corners has all its
// unknown corners in edges[0].
struct OverlapSplit {
    CornerRange inner;
    CornerRange edges[4];  // the top, bottom, left and right strips
};

inline OverlapSplit overlap_split(const Block& b) {
    const CornerRange u = unknown_corners(b);
    const CornerRange in{b.wall_top ? u.i_first : u.i_first + 1, b.wall_bottom ? u.i_last : u.i_last - 1,
                         b.wall_left ? u.j_first : u.j_first + 1, b.wall_right ? u.j_last : u.j_last - 1};
    const CornerRange none{1, 0, 1, 0};
    if (in.empty()) return OverlapSplit{none, {u, none, none, none}};
    return OverlapSplit{in,
                        {CornerRange{u.i_first, in.i_first - 1, u.j_first, u.j_last},
                         CornerRange{in.i_last + 1, u.i_last, u.j_first, u.j_last},
                         CornerRange{in.i_first, in.i_last, u.j_first, in.j_first - 1},
                         CornerRange{in.i_first, in.i_last, in.j_last + 1, u.j_last}}};
}

class PoissonCG {
public:
    // tol: stop when ||r|| <= tol ||zeta||. max_iterations: give up (throw) after this many.
    PoissonCG(const Block& b, double dx, double dy, double tol, int max_iterations)
        : unknown_(unknown_corners(b)), split_(overlap_split(b)), cx_(1 / (dx * dx)), cy_(1 / (dy * dy)),
          diagonal_(2 / (dx * dx) + 2 / (dy * dy)), tol_(tol), max_iterations_(max_iterations),
          r_(b.ny, b.nx), p_(b.ny, b.nx), ap_(b.ny, b.nx) {}

    // Solve A psi = zeta. On entry psi holds the initial guess (the previous
    // step's psi: a warm start); on return, the solution. The wall corners and
    // wall ghosts of psi must be 0. Throws std::runtime_error when the
    // tolerance is not reached within max_iterations.
    CgStats solve(const Field& zeta, Field& psi, const PoissonComm& comm) {
        const CornerRange u = unknown_;
        comm.exchange(psi);

        // r = zeta - A psi, p = r; and the squared norms of r and zeta.
        double sums[2] = {0.0, 0.0};
        double rr0 = 0.0, bb0 = 0.0;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for reduction(+ : rr0, bb0) schedule(static)
#endif
        for (int i = u.i_first; i <= u.i_last; ++i) {
            const double* z = zeta.row(i);
            const double* xn = psi.row(i - 1);
            const double* x = psi.row(i);
            const double* xs = psi.row(i + 1);
            double* r = r_.row(i);
            double* p = p_.row(i);
            for (int j = u.j_first; j <= u.j_last; ++j) {
                r[j] = z[j] - apply(xn, x, xs, j);
                p[j] = r[j];
                rr0 += r[j] * r[j];
                bb0 += z[j] * z[j];
            }
        }
        sums[0] = rr0;
        sums[1] = bb0;
        comm.sum(sums, 2);
        double rr = sums[0];
        const double bb = sums[1];

        // zeta = 0 everywhere: the solution is psi = 0.
        if (bb == 0.0) {
            set_unknown(psi, 0.0);
            return CgStats{0, 0.0};
        }
        const double target = tol_ * tol_ * bb;
        if (rr <= target) return CgStats{0, std::sqrt(rr / bb)};

        for (int k = 1; k <= max_iterations_; ++k) {
            // Ap = A p and p . Ap. With overlap, A is applied to the inner
            // corners while the halo messages of p travel, and to the edge
            // corners once they have arrived.
            double pap = 0.0;
            if (comm.overlap()) {
                comm.start(p_);
                pap = apply_a(split_.inner);
                comm.finish(p_);
                for (const CornerRange& edge : split_.edges) pap += apply_a(edge);
            } else {
                comm.exchange(p_);
                pap = apply_a(u);
            }
            comm.sum(&pap, 1);
            const double alpha = rr / pap;

            // psi += alpha p, r -= alpha Ap, and r . r
            double rr_new = 0.0;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for reduction(+ : rr_new) schedule(static)
#endif
            for (int i = u.i_first; i <= u.i_last; ++i) {
                const double* p = p_.row(i);
                const double* ap = ap_.row(i);
                double* x = psi.row(i);
                double* r = r_.row(i);
                for (int j = u.j_first; j <= u.j_last; ++j) {
                    x[j] += alpha * p[j];
                    r[j] -= alpha * ap[j];
                    rr_new += r[j] * r[j];
                }
            }
            comm.sum(&rr_new, 1);
            if (rr_new <= target) return CgStats{k, std::sqrt(rr_new / bb)};

            // p = r + beta p
            const double beta = rr_new / rr;
            rr = rr_new;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (int i = u.i_first; i <= u.i_last; ++i) {
                const double* r = r_.row(i);
                double* p = p_.row(i);
                for (int j = u.j_first; j <= u.j_last; ++j) p[j] = r[j] + beta * p[j];
            }
        }
        throw std::runtime_error("the Poisson solver (CG) did not reach the tolerance in " +
                                 std::to_string(max_iterations_) + " iterations; relative residual " +
                                 std::to_string(std::sqrt(rr / bb)) + ". Raise --poisson-max-iter or --poisson-tol.");
    }

private:
    // (A f) at column j of a row, given the rows above (n), at (c) and below (s);
    // neighbours on the walls are 0.
    double apply(const double* n, const double* c, const double* s, int j) const {
        return diagonal_ * c[j] - cx_ * (c[j - 1] + c[j + 1]) - cy_ * (n[j] + s[j]);
    }

    // Ap = A p on the corners of `range`; returns the sum of p * Ap over them.
    double apply_a(const CornerRange& range) {
        double pap = 0.0;
        if (range.empty()) return pap;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for reduction(+ : pap) schedule(static)
#endif
        for (int i = range.i_first; i <= range.i_last; ++i) {
            const double* pn = p_.row(i - 1);
            const double* p = p_.row(i);
            const double* ps = p_.row(i + 1);
            double* ap = ap_.row(i);
            for (int j = range.j_first; j <= range.j_last; ++j) {
                ap[j] = apply(pn, p, ps, j);
                pap += p[j] * ap[j];
            }
        }
        return pap;
    }

    void set_unknown(Field& f, double value) const {
        for (int i = unknown_.i_first; i <= unknown_.i_last; ++i) {
            for (int j = unknown_.j_first; j <= unknown_.j_last; ++j) f(i, j) = value;
        }
    }

    CornerRange unknown_;
    OverlapSplit split_;
    double cx_;
    double cy_;
    double diagonal_;
    double tol_;
    int max_iterations_;
    Field r_;   // residual
    Field p_;   // search direction; needs its halos for A p
    Field ap_;  // A p
};
