#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

using Clock = std::chrono::steady_clock;
using Field = std::unique_ptr<double[]>;
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Options {
    int nx = 32, ny = 32, nz = 32;
    int px = 0, py = 0, max_iter = 5000;
    double lx = 1.0, ly = 1.0, lz = 1.0, tol = 1e-8;
    std::string rhs = "benchmark", exchange = "blocking", solution;
};

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string name = argv[i];
        if (i + 1 == argc) throw std::runtime_error("missing value for " + name);
        const std::string v = argv[++i];
        if (name == "--nx") o.nx = std::stoi(v);
        else if (name == "--ny") o.ny = std::stoi(v);
        else if (name == "--nz") o.nz = std::stoi(v);
        else if (name == "--px") o.px = std::stoi(v);
        else if (name == "--py") o.py = std::stoi(v);
        else if (name == "--lx") o.lx = std::stod(v);
        else if (name == "--ly") o.ly = std::stod(v);
        else if (name == "--lz") o.lz = std::stod(v);
        else if (name == "--tol") o.tol = std::stod(v);
        else if (name == "--max-iter") o.max_iter = std::stoi(v);
        else if (name == "--rhs") o.rhs = v;
        else if (name == "--exchange") o.exchange = v;
        else if (name == "--solution") o.solution = v;
        else throw std::runtime_error("unknown option: " + name);
    }
    const auto positive = [](double x) { return std::isfinite(x) && x > 0; };
    if (o.nx < 1 || o.ny < 1 || o.nz < 1 || o.px < 0 || o.py < 0 ||
        o.max_iter < 1 || !positive(o.lx) || !positive(o.ly) || !positive(o.lz) ||
        !positive(o.tol) || o.tol >= 1 ||
        (o.rhs != "benchmark" && o.rhs != "discrete" && o.rhs != "continuous") ||
        (o.exchange != "blocking" && o.exchange != "overlap"))
        throw std::runtime_error("invalid grid, tolerance, iteration limit, or RHS");
    return o;
}

struct Layout {
    int rank = 0, ranks = 1, nodes = 1, px = 1, py = 1, cx = 0, cy = 0;
    int local_rank = 0, local_ranks = 1;
#ifdef USE_MPI
    MPI_Comm cart = MPI_COMM_NULL;
    int neighbor[4]{};
    MPI_Datatype xface = MPI_DATATYPE_NULL, yface = MPI_DATATYPE_NULL;
#endif

    explicit Layout(const Options& o) {
#ifdef USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);
        int dims[2] = {o.py, o.px};
        if ((dims[0] && (dims[0] > ranks || ranks % dims[0])) ||
            (dims[1] && (dims[1] > ranks || ranks % dims[1])) ||
            (dims[0] && dims[1] && 1LL * dims[0] * dims[1] != ranks))
            throw std::runtime_error("process grid does not match MPI rank count");
        MPI_Dims_create(ranks, 2, dims);
        py = dims[0]; px = dims[1];
#else
        if (o.px > 1 || o.py > 1) throw std::runtime_error("serial process grid must be 1 x 1");
#endif
        if (px > o.nx || py > o.ny) throw std::runtime_error("each rank must own at least one point in x and y");
#ifdef USE_MPI
        int periods[2] = {0, 0}, coords[2];
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);
        MPI_Comm_rank(cart, &rank);
        MPI_Cart_coords(cart, rank, 2, coords);
        cy = coords[0]; cx = coords[1];
        MPI_Cart_shift(cart, 1, 1, &neighbor[0], &neighbor[1]);
        MPI_Cart_shift(cart, 0, 1, &neighbor[2], &neighbor[3]);
#endif
    }

    double sum(double local) const {
#ifdef USE_MPI
        if (ranks > 1) {
            double result;
            MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_SUM, cart);
            return result;
        }
#endif
        return local;
    }

    double maximum(double local) const {
#ifdef USE_MPI
        if (ranks > 1) {
            double result;
            MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX, cart);
            return result;
        }
#endif
        return local;
    }

    double minimum(double local) const {
#ifdef USE_MPI
        if (ranks > 1) {
            double result;
            MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_MIN, cart);
            return result;
        }
#endif
        return local;
    }

    void barrier() const {
#ifdef USE_MPI
        MPI_Barrier(cart);
#endif
    }

    // Vanda has one logical CPU per physical core. Check each observed thread
    // has a one-CPU affinity mask and that there are no duplicates on its node.
    std::string placement(int threads) {
        std::vector<int> cores(threads, -1);
#ifdef __linux__
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            cpu_set_t set;
            CPU_ZERO(&set);
            if (sched_getaffinity(0, sizeof(set), &set) == 0 && CPU_COUNT(&set) == 1)
                cores[tid] = sched_getcpu();
        }
#endif
#ifdef USE_MPI
        MPI_Comm shared;
        MPI_Comm_split_type(cart, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &shared);
        MPI_Comm_rank(shared, &local_rank);
        MPI_Comm_size(shared, &local_ranks);
        nodes = static_cast<int>(sum(local_rank == 0 ? 1.0 : 0.0));
        std::vector<int> counts(local_ranks), offsets(local_ranks);
        MPI_Allgather(&threads, 1, MPI_INT, counts.data(), 1, MPI_INT, shared);
        int total = 0;
        for (int i = 0; i < local_ranks; ++i) { offsets[i] = total; total += counts[i]; }
        std::vector<int> local_cores(total);
        MPI_Allgatherv(cores.data(), threads, MPI_INT, local_cores.data(), counts.data(),
                       offsets.data(), MPI_INT, shared);
        MPI_Comm_free(&shared);
        cores.swap(local_cores);
#endif
#ifndef __linux__
        return "unknown";
#else
        std::sort(cores.begin(), cores.end());
        const bool valid = cores.front() >= 0 && std::adjacent_find(cores.begin(), cores.end()) == cores.end();
        return maximum(valid ? 0.0 : 1.0) == 0 ? "cores" : "invalid";
#endif
    }

    ~Layout() {
#ifdef USE_MPI
        if (xface != MPI_DATATYPE_NULL) MPI_Type_free(&xface);
        if (yface != MPI_DATATYPE_NULL) MPI_Type_free(&yface);
        if (cart != MPI_COMM_NULL) MPI_Comm_free(&cart);
#endif
    }
};

struct Grid {
    Layout& layout;
    int nx, ny, nz, x0, y0;
    std::size_t pitch, plane, size;
    double hx, hy, hz, ax, ay, az, diagonal;

    Grid(const Options& o, Layout& l)
        : layout(l), nx(o.nx / l.px + (l.cx < o.nx % l.px)),
          ny(o.ny / l.py + (l.cy < o.ny % l.py)), nz(o.nz),
          x0(l.cx * (o.nx / l.px) + std::min(l.cx, o.nx % l.px)),
          y0(l.cy * (o.ny / l.py) + std::min(l.cy, o.ny % l.py)),
          pitch(static_cast<std::size_t>(nx) + 2), plane(pitch * (static_cast<std::size_t>(ny) + 2)),
          size(plane * (static_cast<std::size_t>(nz) + 2)),
          hx(o.lx / (o.nx + 1.0)), hy(o.ly / (o.ny + 1.0)), hz(o.lz / (o.nz + 1.0)),
          ax(1 / (hx * hx)), ay(1 / (hy * hy)), az(1 / (hz * hz)), diagonal(2 * (ax + ay + az)) {
#ifdef USE_MPI
        if (l.ranks == 1) return;
        MPI_Datatype row_slice;
        MPI_Type_vector(ny, 1, static_cast<int>(pitch), MPI_DOUBLE, &row_slice);
        const auto stride = static_cast<MPI_Aint>(plane * sizeof(double));
        MPI_Type_create_hvector(nz, 1, stride, row_slice, &l.xface);
        MPI_Type_create_hvector(nz, nx, stride, MPI_DOUBLE, &l.yface);
        MPI_Type_commit(&l.xface);
        MPI_Type_commit(&l.yface);
        MPI_Type_free(&row_slice);
#endif
    }

    std::size_t index(int i, int j, int k) const { return i + j * pitch + k * plane; }

#ifdef USE_MPI
    void start_exchange(Field& u, MPI_Request (&requests)[8]) const {
        auto buf = u.get();
        const std::size_t owned[] = {index(1, 1, 1), index(nx, 1, 1), index(1, 1, 1), index(1, ny, 1)};
        const std::size_t ghost[] = {index(0, 1, 1), index(nx+1, 1, 1), index(1, 0, 1), index(1, ny+1, 1)};
        for (int d = 0; d < 4; ++d) {
            const auto type = d < 2 ? layout.xface : layout.yface;
            MPI_Irecv(buf + ghost[d], 1, type, layout.neighbor[d], d ^ 1, layout.cart, &requests[d]);
            MPI_Isend(buf + owned[d], 1, type, layout.neighbor[d], d, layout.cart, &requests[d + 4]);
        }
    }
#endif

    void stencil(const Field& u, Field& out, int i0, int i1, int j0, int j1) const {
        if (i0 > i1 || j0 > j1) return;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int k = 1; k <= nz; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    const auto q = index(i, j, k);
                    out[q] = diagonal * u[q] - ax * (u[q-1] + u[q+1])
                             - ay * (u[q-pitch] + u[q+pitch]) - az * (u[q-plane] + u[q+plane]);
                }
    }

    void apply(Field& u, Field& out, double& halo_time, bool overlap = false) const {
#ifdef USE_MPI
        if (layout.ranks > 1) {
            MPI_Request requests[8];
            const auto start = Clock::now();
            start_exchange(u, requests);
            if (overlap) {
                halo_time += seconds(start);
                const int i0 = layout.neighbor[0] == MPI_PROC_NULL ? 1 : 2;
                const int i1 = layout.neighbor[1] == MPI_PROC_NULL ? nx : nx - 1;
                const int j0 = layout.neighbor[2] == MPI_PROC_NULL ? 1 : 2;
                const int j1 = layout.neighbor[3] == MPI_PROC_NULL ? ny : ny - 1;
                const bool has_interior = i0 <= i1 && j0 <= j1;
                if (has_interior) stencil(u, out, i0, i1, j0, j1);
                const auto wait_start = Clock::now();
                MPI_Waitall(8, requests, MPI_STATUSES_IGNORE);
                halo_time += seconds(wait_start);
                if (has_interior) {
                    stencil(u, out, 1, nx, 1, j0 - 1);
                    stencil(u, out, 1, nx, j1 + 1, ny);
                    stencil(u, out, 1, i0 - 1, j0, j1);
                    stencil(u, out, i1 + 1, nx, j0, j1);
                    return;
                }
            } else {
                MPI_Waitall(8, requests, MPI_STATUSES_IGNORE);
                halo_time += seconds(start);
            }
        }
#else
        (void)halo_time;
        (void)overlap;
#endif
        stencil(u, out, 1, nx, 1, ny);
    }

    double dot(const Field& a, const Field& b) const {
        double local = 0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) reduction(+:local)
#endif
        for (int k = 1; k <= nz; ++k)
            for (int j = 1; j <= ny; ++j)
                for (int i = 1; i <= nx; ++i) {
                    const auto q = index(i, j, k);
                    local += a[q] * b[q];
                }
        return local;
    }
};

double pattern(std::uint64_t key) {
    key += 0x9e3779b97f4a7c15ULL;
    key = (key ^ (key >> 30)) * 0xbf58476d1ce4e5b9ULL;
    key = (key ^ (key >> 27)) * 0x94d049bb133111ebULL;
    return static_cast<double>((key ^ (key >> 31)) >> 11) * 0x1.0p-53 - 0.5;
}

struct Problem {
    const Options& o;
    const Grid& g;
    const double pi = std::acos(-1);
    const int modes[2][3] = {{1, 1, 1}, {2, 3, 4}};
    double lambda[2]{};

    Problem(const Options& opts, const Grid& grid) : o(opts), g(grid) {
        const double lengths[] = {o.lx, o.ly, o.lz}, invh2[] = {g.ax, g.ay, g.az};
        const int n[] = {o.nx, o.ny, o.nz};
        for (int m = 0; m < 2; ++m)
            for (int d = 0; d < 3; ++d) {
                if (o.rhs == "discrete") {
                    const double s = std::sin(pi * modes[m][d] / (2.0 * (n[d] + 1.0)));
                    lambda[m] += 4 * invh2[d] * s * s;
                } else {
                    const double wave = pi * modes[m][d] / lengths[d];
                    lambda[m] += wave * wave;
                }
            }
    }

    void terms(int i, int j, int k, double& a, double& b) const {
        const double x = pi * (g.x0 + i) / (o.nx + 1.0);
        const double y = pi * (g.y0 + j) / (o.ny + 1.0);
        const double z = pi * k / (o.nz + 1.0);
        a = std::sin(x) * std::sin(y) * std::sin(z);
        b = 0.1 * std::sin(2*x) * std::sin(3*y) * std::sin(4*z);
    }

    std::size_t global_index(int i, int j, int k) const {
        return (static_cast<std::size_t>(k-1)*o.ny + g.y0+j-1)*o.nx + g.x0+i-1;
    }

    double rhs(int i, int j, int k) const {
        double a, b;
        terms(i, j, k, a, b);
        if (o.rhs != "benchmark") return lambda[0]*a + lambda[1]*b;
        // Hash the global coordinates, so rank layout does not change the RHS.
        const auto key = (std::uint64_t(g.x0+i) << 42) ^ (std::uint64_t(g.y0+j) << 21) ^ std::uint64_t(k);
        return a + pattern(key);
    }
};

void save_solution(const Options& o, const Grid& g, const Problem& problem,
                   const Field& u, std::size_t n) {
    if (o.solution.empty()) return;
    std::vector<double> values(n, 0.0);
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i)
                values[problem.global_index(i, j, k)] = u[g.index(i, j, k)];
#ifdef USE_MPI
    if (g.layout.rank == 0)
        MPI_Reduce(MPI_IN_PLACE, values.data(), static_cast<int>(n), MPI_DOUBLE, MPI_SUM, 0, g.layout.cart);
    else
        MPI_Reduce(values.data(), nullptr, static_cast<int>(n), MPI_DOUBLE, MPI_SUM, 0, g.layout.cart);
#endif
    if (g.layout.rank == 0) {
        std::ofstream out(o.solution);
        out << std::setprecision(17);
        for (double v : values) out << v << '\n';
        if (!out) throw std::runtime_error("cannot write " + o.solution);
    }
}

#ifdef USE_GPU
#include "poisson3d_gpu.hpp"
#endif

int run(const Options& o) {
#ifdef USE_GPU
    if (o.exchange != "blocking") throw std::runtime_error("GPU exchange must be blocking");
#endif
    int threads = 1;
#ifdef _OPENMP
    omp_set_dynamic(0);
#pragma omp parallel
    {
#pragma omp single
        threads = omp_get_num_threads();
    }
#endif
    const std::size_t n = static_cast<std::size_t>(o.nx) * o.ny * o.nz;
    if (!o.solution.empty() && n > 1000000)
        throw std::runtime_error("--solution is for validation grids of up to 1,000,000 unknowns");
    const auto nnz = 7*n - 2*(static_cast<std::size_t>(o.nx)*o.ny
                             + static_cast<std::size_t>(o.nx)*o.nz + static_cast<std::size_t>(o.ny)*o.nz);
    Layout layout(o);
    const auto binding = layout.placement(threads);
#ifdef USE_GPU
    require_gpu_device(layout);
#endif
    const Grid g(o, layout);
    const Problem problem(o, g);

    // Allocate without initialization, then first-touch each row on the
    // thread that will use it. Ghost values at physical walls remain zero.
    const auto alloc = [&] { return Field(new double[g.size]); };
    Field b = alloc(), r = alloc(), u = alloc(), p = alloc(), ap = alloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int k = 0; k <= g.nz + 1; ++k)
        for (int j = 0; j <= g.ny + 1; ++j)
            for (int i = 0; i <= g.nx + 1; ++i) {
                const auto q = g.index(i, j, k);
                b[q] = r[q] = u[q] = p[q] = ap[q] = 0;
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                const auto q = g.index(i, j, k);
                b[q] = r[q] = p[q] = problem.rhs(i, j, k);
            }

    double halo_s = 0, reduce_s = 0;
    double device_setup_s = 0, device_finish_s = 0, halo_device_copy_s = 0;
    int gpu_count = 0;
#ifdef USE_GPU
    layout.barrier();
    const auto gpu = solve_gpu(o, g, r, u, p, ap);
    const double bb = gpu.bb, solve_s = layout.maximum(gpu.solve_s);
    const double solve_total_s = layout.maximum(gpu.total_s);
    const int iterations = gpu.iterations;
    device_setup_s = layout.maximum(gpu.setup_s);
    device_finish_s = layout.maximum(gpu.finish_s);
    halo_device_copy_s = layout.maximum(gpu.halo_copy_s);
    halo_s = gpu.halo_s;
    reduce_s = gpu.reduction_s;
    gpu_count = layout.ranks;
    const char* poisson_device = "gpu";
#else
    auto global_dot = [&](double local) {
        if (layout.ranks == 1) return local;
        const auto start = Clock::now();
        const auto global = layout.sum(local);
        reduce_s += seconds(start);
        return global;
    };

    layout.barrier();
    const auto start = Clock::now();
    const double bb = global_dot(g.dot(b, b));
    if (!(bb > 0) || !std::isfinite(bb)) throw std::runtime_error("invalid or zero RHS norm");
    double rr = bb;
    const double target = o.tol * o.tol * bb;
    int iterations = 0;
    while (iterations < o.max_iter) {
        g.apply(p, ap, halo_s, o.exchange == "overlap");
        const double pap = global_dot(g.dot(p, ap));
        if (!(pap > 0) || !std::isfinite(pap)) throw std::runtime_error("CG breakdown");
        const double alpha = rr / pap;
        double local_rr = 0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) reduction(+:local_rr)
#endif
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    const auto q = g.index(i, j, k);
                    u[q] += alpha * p[q];
                    r[q] -= alpha * ap[q];
                    local_rr += r[q] * r[q];
                }
        const double next_rr = global_dot(local_rr);
        ++iterations;
        if (!std::isfinite(next_rr)) throw std::runtime_error("CG residual is not finite");
        if (next_rr <= target) break;
        const double beta = next_rr / rr;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
        for (int k = 1; k <= g.nz; ++k)
            for (int j = 1; j <= g.ny; ++j)
                for (int i = 1; i <= g.nx; ++i) {
                    const auto q = g.index(i, j, k);
                    p[q] = r[q] + beta * p[q];
                }
        rr = next_rr;
    }
    const double solve_s = layout.maximum(seconds(start));
    const double solve_total_s = solve_s;
    const char* poisson_device = "cpu";
#endif
    halo_s = layout.maximum(halo_s);
    const double reduction_min_s = layout.minimum(reduce_s);
    reduce_s = layout.maximum(reduce_s);

    double ignored = 0;
    g.apply(u, ap, ignored);
    double residual2 = 0, error2 = 0, exact2 = 0, max_error = 0;
    for (int k = 1; k <= g.nz; ++k)
        for (int j = 1; j <= g.ny; ++j)
            for (int i = 1; i <= g.nx; ++i) {
                const auto q = g.index(i, j, k);
                const double residual = b[q] - ap[q];
                residual2 += residual * residual;
                if (o.rhs != "benchmark") {
                    double a, term;
                    problem.terms(i, j, k, a, term);
                    const double exact = a + term, delta = u[q] - exact;
                    error2 += delta * delta;
                    exact2 += exact * exact;
                    max_error = std::max(max_error, std::abs(delta));
                }
            }
    const double rel_residual = std::sqrt(layout.sum(residual2) / bb);
    if (!(rel_residual <= 1.1 * o.tol)) throw std::runtime_error("CG did not attain the requested true-residual tolerance");
    const double rel_error = o.rhs == "benchmark" ? 0.0 : std::sqrt(layout.sum(error2) / layout.sum(exact2));
    max_error = layout.maximum(max_error);
    save_solution(o, g, problem, u, n);

    if (layout.rank == 0) {
        std::cout << std::setprecision(17)
                  << "{\"nx\":" << o.nx << ",\"ny\":" << o.ny << ",\"nz\":" << o.nz
                  << ",\"lx\":" << o.lx << ",\"ly\":" << o.ly << ",\"lz\":" << o.lz
                  << ",\"px\":" << layout.px << ",\"py\":" << layout.py
                  << ",\"ranks\":" << layout.ranks << ",\"threads\":" << threads
                  << ",\"nodes\":" << layout.nodes << ",\"binding\":\"" << binding
                  << "\",\"tol\":" << o.tol << ",\"rhs\":\"" << o.rhs << "\",\"N\":" << n
                  << ",\"nnz\":" << nnz << ",\"iterations\":" << iterations
                  << ",\"solve_s\":" << solve_s << ",\"halo_s\":" << halo_s
                  << ",\"reduction_s\":" << reduce_s << ",\"reduction_min_s\":" << reduction_min_s
                  << ",\"relative_residual\":" << rel_residual
                  << ",\"poisson_device\":\"" << poisson_device << "\""
                  << ",\"gpu_count\":" << gpu_count << ",\"exchange\":\"" << o.exchange << "\""
                  << ",\"device_setup_s\":" << device_setup_s << ",\"device_finish_s\":" << device_finish_s
                  << ",\"halo_device_copy_s\":" << halo_device_copy_s
                  << ",\"solve_total_s\":" << solve_total_s;
        if (o.rhs == "benchmark") std::cout << ",\"relative_error\":null,\"max_error\":null";
        else std::cout << ",\"relative_error\":" << rel_error << ",\"max_error\":" << max_error;
        std::cout << "}\n";
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "poisson3d [--nx N --ny N --nz N] [--lx L --ly L --lz L]\n"
                     "          [--px P --py P] [--rhs benchmark|discrete|continuous]\n"
                     "          [--exchange blocking|overlap] [--tol EPS --max-iter N] [--solution FILE]\n";
        return 0;
    }
#ifdef USE_MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    int status = 1;
    try {
        status = run(parse(argc, argv));
    } catch (const std::exception& e) {
        std::cerr << "poisson3d: " << e.what() << '\n';
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return status;
}
