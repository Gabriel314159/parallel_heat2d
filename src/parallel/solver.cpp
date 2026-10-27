#include "parallel/solver.hpp"
#include "parallel/halo.hpp"
#include "common/stencil.hpp"
#include <cstdlib>
#ifdef HEAT2D_USE_OPENMP
#include <omp.h>
#endif

namespace heat2d {
int configure_threads(int requested) {
    int actual=1;
#ifdef HEAT2D_USE_OPENMP
    omp_set_dynamic(0);
    if (requested>0) omp_set_num_threads(requested);
    else if (!std::getenv("OMP_NUM_THREADS")) omp_set_num_threads(1);
#pragma omp parallel
    {
#pragma omp single
        actual=omp_get_num_threads();
    }
#else
    if (requested>1) throw std::invalid_argument("This executable was built without OpenMP.");
#endif
    return actual;
}
namespace {
void update_full(const Grid& old, Grid& next, const Config& c) {
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i=1; i<=old.ny; ++i)
        for (int j=1; j<=old.nx; ++j) next(i,j)=update_point(old,i,j,c.rx,c.ry);
}
void update_interior(const Grid& old, Grid& next, const Config& c) {
    if (old.nx<3 || old.ny<3) return;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i=2; i<old.ny; ++i)
        for (int j=2; j<old.nx; ++j) next(i,j)=update_point(old,i,j,c.rx,c.ry);
}
void update_edges(const Grid& old, Grid& next, const Config& c) {
    // Top/bottom rows and the remaining side columns are disjoint, even for 1-cell blocks.
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel
#endif
    {
#ifdef HEAT2D_USE_OPENMP
#pragma omp for schedule(static) nowait
#endif
        for (int j=1; j<=old.nx; ++j) {
            next(1,j)=update_point(old,1,j,c.rx,c.ry);
            if (old.ny>1) next(old.ny,j)=update_point(old,old.ny,j,c.rx,c.ry);
        }
#ifdef HEAT2D_USE_OPENMP
#pragma omp for schedule(static)
#endif
        for (int i=2; i<old.ny; ++i) {
            next(i,1)=update_point(old,i,1,c.rx,c.ry);
            if (old.nx>1) next(i,old.nx)=update_point(old,i,old.nx,c.rx,c.ry);
        }
    }
}
}
double solve_parallel(Grid& current, Grid& next, const Config& c, const Domain& d) {
    Halo halo(d);
    mpi_check(MPI_Barrier(d.comm),"MPI_Barrier before timing");
    const double start=MPI_Wtime();
    for (int step=0; step<c.steps; ++step) {
        if (c.exchange=="blocking") {
            halo.blocking(current);
            update_full(current,next,c);
        } else {
            halo.start(current);
            update_interior(current,next,c);
            halo.finish();
            update_edges(current,next,c);
        }
        current.swap_data(next);
    }
    const double local=MPI_Wtime()-start;
    double maximum=0;
    mpi_check(MPI_Reduce(&local,&maximum,1,MPI_DOUBLE,MPI_MAX,0,d.comm),"MPI_Reduce elapsed");
    return maximum;
}
}
