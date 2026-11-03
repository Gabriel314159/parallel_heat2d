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
//
// With HEAT2D_USE_OFFLOAD the solve runs on a GPU (docs/design.md, "The
// Poisson solve on a GPU").

#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef HEAT2D_USE_OFFLOAD
#include <omp.h>
#endif

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
    // False if no side of the block borders another rank.
    virtual bool has_neighbours() const { return true; }
};

// One process: local sums are global sums, and there are no process halos.
class SerialComm : public PoissonComm {
public:
    void sum(double*, int) const override {}
    void exchange(Field&) const override {}
    bool has_neighbours() const override { return false; }
};

struct CgStats {
    int iterations = 0;
    double relative_residual = 0.0;  // ||r|| / ||zeta|| at the end
    double copy_seconds = 0.0;       // GPU build: host-GPU copies
    bool on_device = false;          // GPU build: the loops ran on a GPU
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

#ifdef HEAT2D_USE_OFFLOAD
// Index of entry q of the edge buffers: rows 1, ny, then columns 1, nx; with
// `ghost`, the ghost next to it (rows 0, ny + 1, columns 0, nx + 1).
#pragma omp declare target
inline std::size_t edge_index(int q, int s, int nx, int ny, bool ghost) {
    const int out = ghost ? 1 : 0;
    int i, j;
    if (q < nx) {
        i = 1 - out;
        j = q + 1;
    } else if (q < 2 * nx) {
        i = ny + out;
        j = q - nx + 1;
    } else if (q < 2 * nx + ny) {
        i = q - 2 * nx + 1;
        j = 1 - out;
    } else {
        i = q - 2 * nx - ny + 1;
        j = nx + out;
    }
    return static_cast<std::size_t>(i) * s + j;
}
#pragma omp end declare target
#endif

class PoissonCG {
public:
    // tol: stop when ||r|| <= tol ||zeta||. max_iterations: give up (throw) after this many.
    PoissonCG(const Block& b, double dx, double dy, double tol, int max_iterations)
        : unknown_(unknown_corners(b)), split_(overlap_split(b)), cx_(1 / (dx * dx)), cy_(1 / (dy * dy)),
          diagonal_(2 / (dx * dx) + 2 / (dy * dy)), tol_(tol), max_iterations_(max_iterations),
          r_(b.ny, b.nx), p_(b.ny, b.nx), ap_(b.ny, b.nx),
          edges_(static_cast<std::size_t>(2 * (b.nx + b.ny))), ghosts_(edges_.size()) {
#ifdef HEAT2D_USE_OFFLOAD
        // The work fields stay on the GPU for the solver's lifetime.
        double* r = r_.data();
        double* p = p_.data();
        double* ap = ap_.data();
        double* e = edges_.data();
        double* g = ghosts_.data();
        const std::size_t n = r_.size();
        const std::size_t m = edges_.size();
#pragma omp target enter data map(to : r[0:n], p[0:n], ap[0:n], e[0:m], g[0:m])
#endif
    }

#ifdef HEAT2D_USE_OFFLOAD
    ~PoissonCG() {
        double* r = r_.data();
        double* p = p_.data();
        double* ap = ap_.data();
        double* e = edges_.data();
        double* g = ghosts_.data();
        const std::size_t n = r_.size();
        const std::size_t m = edges_.size();
#pragma omp target exit data map(delete : r[0:n], p[0:n], ap[0:n], e[0:m], g[0:m])
    }
    // A copy would share the GPU fields.
    PoissonCG(const PoissonCG&) = delete;
    PoissonCG& operator=(const PoissonCG&) = delete;
#endif

    // Solve A psi = zeta. On entry psi holds the initial guess (the previous
    // step's psi: a warm start); on return, the solution. The wall corners and
    // wall ghosts of psi must be 0 and its process halos current (psi is not
    // exchanged here). Throws std::runtime_error when the tolerance is not
    // reached within max_iterations.
    CgStats solve(const Field& zeta, Field& psi, const PoissonComm& comm) {
#ifdef HEAT2D_USE_OFFLOAD
        return solve_on_device(zeta, psi, comm);
#else
        const CornerRange u = unknown_;

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
#endif
    }

private:
#ifdef HEAT2D_USE_OFFLOAD
    // solve() with its loops on the GPU; sums and exchanges stay on the host.
    // Entry (i, j) is at i * s + j.
    CgStats solve_on_device(const Field& zeta, Field& psi, const PoissonComm& comm) {
        static_assert(halo == 1, "the GPU loops read the neighbours at k +- 1 and k +- s");
        const int i0 = unknown_.i_first, i1 = unknown_.i_last;
        const int j0 = unknown_.j_first, j1 = unknown_.j_last;
        const int s = r_.stride();
        const std::size_t n = r_.size();
        const double* z = zeta.data();
        double* x = psi.data();
        double* r = r_.data();
        double* p = p_.data();
        double* ap = ap_.data();
        const double cx = cx_, cy = cy_, dg = diagonal_;
        CgStats stats;

        double start = omp_get_wtime();
#pragma omp target enter data map(to : z[0:n], x[0:n])
        stats.copy_seconds += omp_get_wtime() - start;
        int on_host = 1;
#pragma omp target map(from : on_host)
        on_host = omp_is_initial_device();
        stats.on_device = on_host == 0;

        // r = zeta - A psi, p = r; and the squared norms of r and zeta.
        double rr0 = 0.0, bb0 = 0.0;
#pragma omp target teams distribute parallel for collapse(2) reduction(+ : rr0, bb0) map(tofrom : rr0, bb0)
        for (int i = i0; i <= i1; ++i) {
            for (int j = j0; j <= j1; ++j) {
                const std::size_t k = static_cast<std::size_t>(i) * s + j;
                const double rk = z[k] - (dg * x[k] - cx * (x[k - 1] + x[k + 1]) - cy * (x[k - s] + x[k + s]));
                r[k] = rk;
                p[k] = rk;
                rr0 += rk * rk;
                bb0 += z[k] * z[k];
            }
        }
        double sums[2] = {rr0, bb0};
        comm.sum(sums, 2);
        double rr = sums[0];
        const double bb = sums[1];

        // zeta = 0 everywhere: the solution is psi = 0.
        if (bb == 0.0) {
#pragma omp target teams distribute parallel for collapse(2)
            for (int i = i0; i <= i1; ++i) {
                for (int j = j0; j <= j1; ++j) x[static_cast<std::size_t>(i) * s + j] = 0.0;
            }
            release_device(z, x, n, stats);
            return stats;
        }
        const double target = tol_ * tol_ * bb;
        if (rr <= target) {
            release_device(z, x, n, stats);
            stats.relative_residual = std::sqrt(rr / bb);
            return stats;
        }

        for (int it = 1; it <= max_iterations_; ++it) {
            exchange_on_device(comm, stats);
            double pap = 0.0;
#pragma omp target teams distribute parallel for collapse(2) reduction(+ : pap) map(tofrom : pap)
            for (int i = i0; i <= i1; ++i) {
                for (int j = j0; j <= j1; ++j) {
                    const std::size_t k = static_cast<std::size_t>(i) * s + j;
                    ap[k] = dg * p[k] - cx * (p[k - 1] + p[k + 1]) - cy * (p[k - s] + p[k + s]);
                    pap += p[k] * ap[k];
                }
            }
            comm.sum(&pap, 1);
            const double alpha = rr / pap;

            double rr_new = 0.0;
#pragma omp target teams distribute parallel for collapse(2) reduction(+ : rr_new) map(tofrom : rr_new)
            for (int i = i0; i <= i1; ++i) {
                for (int j = j0; j <= j1; ++j) {
                    const std::size_t k = static_cast<std::size_t>(i) * s + j;
                    x[k] += alpha * p[k];
                    r[k] -= alpha * ap[k];
                    rr_new += r[k] * r[k];
                }
            }
            comm.sum(&rr_new, 1);
            if (rr_new <= target) {
                release_device(z, x, n, stats);
                stats.iterations = it;
                stats.relative_residual = std::sqrt(rr_new / bb);
                return stats;
            }

            const double beta = rr_new / rr;
            rr = rr_new;
#pragma omp target teams distribute parallel for collapse(2)
            for (int i = i0; i <= i1; ++i) {
                for (int j = j0; j <= j1; ++j) {
                    const std::size_t k = static_cast<std::size_t>(i) * s + j;
                    p[k] = r[k] + beta * p[k];
                }
            }
        }
        release_device(z, x, n, stats);
        throw std::runtime_error("the Poisson solver (CG) did not reach the tolerance in " +
                                 std::to_string(max_iterations_) + " iterations; relative residual " +
                                 std::to_string(std::sqrt(rr / bb)) + ". Raise --poisson-max-iter or --poisson-tol.");
    }

    // psi back to the host, zeta released, on every way out of a solve. Not a
    // lambda: nvc++ 26.1 skipped this directive inside one.
    static void release_device(const double* z, double* x, std::size_t n, CgStats& stats) {
        const double t = omp_get_wtime();
#pragma omp target exit data map(from : x[0:n]) map(release : z[0:n])
        stats.copy_seconds += omp_get_wtime() - t;
    }

    // Fill the process halos of p on the GPU: only its edges go to the host
    // copy of p, which is exchanged as usual, and only the ghosts come back.
    void exchange_on_device(const PoissonComm& comm, CgStats& stats) {
        if (!comm.has_neighbours()) return;
        const int s = p_.stride();
        const int nx = p_.nx(), ny = p_.ny();
        const int m = static_cast<int>(edges_.size());
        double* p = p_.data();
        double* e = edges_.data();
        double* g = ghosts_.data();
#pragma omp target teams distribute parallel for
        for (int q = 0; q < m; ++q) {
            e[q] = p[edge_index(q, s, nx, ny, false)];
        }
        double start = omp_get_wtime();
#pragma omp target update from(e[0:m])
        stats.copy_seconds += omp_get_wtime() - start;
        for (int q = 0; q < m; ++q) p[edge_index(q, s, nx, ny, false)] = e[q];
        comm.exchange(p_);
        for (int q = 0; q < m; ++q) g[q] = p[edge_index(q, s, nx, ny, true)];
        start = omp_get_wtime();
#pragma omp target update to(g[0:m])
        stats.copy_seconds += omp_get_wtime() - start;
#pragma omp target teams distribute parallel for
        for (int q = 0; q < m; ++q) {
            p[edge_index(q, s, nx, ny, true)] = g[q];
        }
    }

#endif

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
    // GPU build: the edges and ghosts of p in transit (exchange_on_device).
    std::vector<double> edges_;
    std::vector<double> ghosts_;
};
