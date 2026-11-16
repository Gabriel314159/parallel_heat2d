#pragma once

#if !defined(_OPENMP)
#error The GPU solver requires OpenMP target
#endif

struct GpuSolve {
    double bb = 0;
    int iterations = 0;
    double solve_s = 0, setup_s = 0, finish_s = 0, total_s = 0;
    double halo_s = 0, halo_copy_s = 0, reduction_s = 0;
};

void require_gpu_device(const Layout& layout) {
    if (layout.sum(omp_get_num_devices() < layout.local_ranks ? 1.0 : 0.0) != 0)
        throw std::runtime_error("fewer visible OpenMP devices than MPI ranks on a node");
    omp_set_default_device(layout.local_rank);
    int on_host = 1;
#pragma omp target map(from:on_host)
    { on_host = omp_is_initial_device(); }
    if (layout.sum(on_host) != 0) throw std::runtime_error("OpenMP target fell back to the host");
}

#ifdef USE_MPI
#pragma omp declare target
static std::size_t face_index(std::size_t q, std::size_t xf, std::size_t yf,
                              int nx, int ny, std::size_t s, std::size_t plane, bool ghost) {
    // West/east faces contain ny*nz values; south/north contain nx*nz.
    if (q < 2*xf) {
        const bool east = q >= xf;
        q %= xf;
        const int i = east ? nx + ghost : 1 - ghost;
        return (q / ny + 1) * plane + (q % ny + 1) * s + i;
    }
    q -= 2*xf;
    const bool north = q >= yf;
    q %= yf;
    const int j = north ? ny + ghost : 1 - ghost;
    return (q / nx + 1) * plane + j * s + q % nx + 1;
}
#pragma omp end declare target

static void exchange_gpu_direction(double* p, double* send, double* recv, const Grid& g, GpuSolve& stats) {
    const int nx = g.nx, ny = g.ny;
    const auto n = g.size, s = g.pitch, plane = g.plane;
    const auto xf = static_cast<std::size_t>(ny) * g.nz, yf = static_cast<std::size_t>(nx) * g.nz;
    const auto size = 2 * (xf + yf);
    const auto halo_start = Clock::now();
#pragma omp target teams distribute parallel for map(alloc:p[0:n], send[0:size])
    for (std::size_t q = 0; q < size; ++q)
        send[q] = p[face_index(q, xf, yf, nx, ny, s, plane, false)];

    const auto to_host = Clock::now();
#pragma omp target update from(send[0:size])
    stats.halo_copy_s += seconds(to_host);
    const std::size_t offsets[] = {0, xf, 2*xf, 2*xf+yf};
    MPI_Request requests[8];
    for (int d = 0; d < 4; ++d) {
        const int count = static_cast<int>(d < 2 ? xf : yf);
        const auto peer = g.layout.neighbor[d];
        MPI_Irecv(recv + offsets[d], count, MPI_DOUBLE, peer, d ^ 1, g.layout.cart, &requests[d]);
        MPI_Isend(send + offsets[d], count, MPI_DOUBLE, peer, d, g.layout.cart, &requests[d+4]);
    }
    MPI_Waitall(8, requests, MPI_STATUSES_IGNORE);
    // MPI_PROC_NULL leaves the initially zero receive entries untouched.
    const auto to_device = Clock::now();
#pragma omp target update to(recv[0:size])
    stats.halo_copy_s += seconds(to_device);
#pragma omp target teams distribute parallel for map(alloc:p[0:n], recv[0:size])
    for (std::size_t q = 0; q < size; ++q)
        p[face_index(q, xf, yf, nx, ny, s, plane, true)] = recv[q];
    stats.halo_s += seconds(halo_start);
}
#endif

GpuSolve solve_gpu(const Options& o, const Grid& g, Field& residual, Field& solution,
                   Field& direction, Field& product) {
    const auto n = g.size, s = g.pitch, plane = g.plane;
    const int nx = g.nx, ny = g.ny, nz = g.nz;
    const double ax = g.ax, ay = g.ay, az = g.az, dg = g.diagonal;
    double* r = residual.get();
    double* x = solution.get();
    double* p = direction.get();
    double* ap = product.get();

    GpuSolve stats;
#ifdef USE_MPI
    const bool distributed = g.layout.ranks > 1;
    const std::size_t face_size = distributed ? 2 * (static_cast<std::size_t>(nx) + ny) * nz : 0;
    std::vector<double> send_buffer(face_size), recv_buffer(face_size, 0.0);
    double* send = send_buffer.data();
    double* recv = recv_buffer.data();
#endif
    const auto total_start = Clock::now();
    // Physical wall values in p and x are zero. Keep all arrays resident
    // during iteration. Only the reduction scalars come back each round.
#pragma omp target enter data map(to:r[0:n], x[0:n], p[0:n]) map(alloc:ap[0:n])
#ifdef USE_MPI
    if (distributed) {
#pragma omp target enter data map(alloc:send[0:face_size], recv[0:face_size])
    }
#endif
    stats.setup_s = seconds(total_start);
    const auto solve_start = Clock::now();

    double bb = 0;
#pragma omp target teams distribute parallel for collapse(3) map(alloc:r[0:n]) reduction(+:bb) map(tofrom:bb)
    for (int k = 1; k <= nz; ++k)
        for (int j = 1; j <= ny; ++j)
            for (int i = 1; i <= nx; ++i) {
                const std::size_t q = k * plane + j * s + i;
                bb += r[q] * r[q];
            }
    const auto global_dot = [&](double local) {
        if (g.layout.ranks == 1) return local;
        const auto start = Clock::now();
        const auto result = g.layout.sum(local);
        stats.reduction_s += seconds(start);
        return result;
    };
    bb = global_dot(bb);
    stats.bb = bb;
    const char* failure = nullptr;
    if (!(bb > 0) || !std::isfinite(bb)) failure = "invalid or zero RHS norm";

    double rr = bb;
    const double target = o.tol * o.tol * bb;
    while (!failure && stats.iterations < o.max_iter) {
#ifdef USE_MPI
        if (distributed) exchange_gpu_direction(p, send, recv, g, stats);
#endif
        double pap = 0;
#pragma omp target teams distribute parallel for collapse(3) map(alloc:p[0:n], ap[0:n]) reduction(+:pap) map(tofrom:pap)
        for (int k = 1; k <= nz; ++k)
            for (int j = 1; j <= ny; ++j)
                for (int i = 1; i <= nx; ++i) {
                    const std::size_t q = k * plane + j * s + i;
                    const double value = dg * p[q] - ax * (p[q-1] + p[q+1])
                                         - ay * (p[q-s] + p[q+s]) - az * (p[q-plane] + p[q+plane]);
                    ap[q] = value;
                    pap += p[q] * value;
                }
        pap = global_dot(pap);
        if (!(pap > 0) || !std::isfinite(pap)) {
            failure = "GPU CG breakdown";
            break;
        }

        const double alpha = rr / pap;
        double next_rr = 0;
#pragma omp target teams distribute parallel for collapse(3) map(alloc:x[0:n], r[0:n], p[0:n], ap[0:n]) reduction(+:next_rr) map(tofrom:next_rr)
        for (int k = 1; k <= nz; ++k)
            for (int j = 1; j <= ny; ++j)
                for (int i = 1; i <= nx; ++i) {
                    const std::size_t q = k * plane + j * s + i;
                    x[q] += alpha * p[q];
                    r[q] -= alpha * ap[q];
                    next_rr += r[q] * r[q];
                }
        next_rr = global_dot(next_rr);
        ++stats.iterations;
        if (!std::isfinite(next_rr)) {
            failure = "GPU CG residual is not finite";
            break;
        }
        if (next_rr <= target) break;

        const double beta = next_rr / rr;
#pragma omp target teams distribute parallel for collapse(3) map(alloc:p[0:n], r[0:n])
        for (int k = 1; k <= nz; ++k)
            for (int j = 1; j <= ny; ++j)
                for (int i = 1; i <= nx; ++i) {
                    const std::size_t q = k * plane + j * s + i;
                    p[q] = r[q] + beta * p[q];
                }
        rr = next_rr;
    }
    stats.solve_s = seconds(solve_start);

    const auto finish_start = Clock::now();
#pragma omp target exit data map(from:x[0:n]) map(release:r[0:n], p[0:n], ap[0:n])
#ifdef USE_MPI
    if (distributed) {
#pragma omp target exit data map(release:send[0:face_size], recv[0:face_size])
    }
#endif
    stats.finish_s = seconds(finish_start);
    stats.total_s = seconds(total_start);
    if (failure) throw std::runtime_error(failure);
    return stats;
}
