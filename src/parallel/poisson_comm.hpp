// The global sums and the halo exchanges of the CG Poisson solver over MPI,
// and the time spent in them. See docs/design.md, "Poisson equation for ψ"
// and "Output and metrics".

#pragma once

#include <mpi.h>

#include "common/grid.hpp"
#include "common/poisson.hpp"
#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"

// Seconds spent inside the communication calls of the CG solver.
struct CommTimes {
    // Halo exchanges: the whole exchange when blocking; when overlapping,
    // only posting the messages and waiting for them, the part that the
    // computation did not hide.
    double halo = 0.0;
    // Global sums, including the wait for the slowest rank.
    double sum = 0.0;
};

class MpiComm : public PoissonComm {
public:
    MpiComm(MPI_Comm comm, HaloExchange& halos, bool overlap) : comm_(comm), halos_(halos), overlap_(overlap) {}

    void sum(double* values, int count) const override {
        const double start = MPI_Wtime();
        mpi_check(MPI_Allreduce(MPI_IN_PLACE, values, count, MPI_DOUBLE, MPI_SUM, comm_), "Poisson solver sum");
        times_.sum += MPI_Wtime() - start;
    }
    void exchange(Field& f) const override {
        const double start = MPI_Wtime();
        halos_.exchange(f);
        times_.halo += MPI_Wtime() - start;
    }
    bool overlap() const override { return overlap_; }
    void start(Field& f) const override {
        const double start = MPI_Wtime();
        halos_.start(f);
        times_.halo += MPI_Wtime() - start;
    }
    void finish(Field& f) const override {
        const double start = MPI_Wtime();
        halos_.finish(f);
        times_.halo += MPI_Wtime() - start;
    }

    // The times since the last reset.
    CommTimes times() const { return times_; }
    void reset_times() { times_ = CommTimes{}; }

private:
    MPI_Comm comm_;
    HaloExchange& halos_;
    bool overlap_;
    // The functions above are const: they do not change what the solver
    // computes. "mutable" lets them add to the times all the same.
    mutable CommTimes times_;
};
