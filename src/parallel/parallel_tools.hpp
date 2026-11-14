// Small MPI and OpenMP helpers shared by the 2D and 3D parallel solvers
// (heat_parallel and heat3d_parallel): global sums and maxima, thread setup,
// the binding warning, rank-0 file output and the strict interior for overlap.
//
// Header-only: each solver compiles it with its own OpenMP setting.

#pragma once

#include <mpi.h>

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef HEAT2D_USE_OPENMP
#include <omp.h>
#endif

#include "common/affinity.hpp"
#include "common/grid.hpp"
#include "parallel/decomposition.hpp"

// Set the number of OpenMP threads of this rank: `requested` (--threads) if it
// is positive, else OMP_NUM_THREADS if that is set, else 1. Dynamic adjustment
// is switched off, so the number is what was asked for. Returns the number of
// threads in use. Throws if more than one thread is requested in a build
// without OpenMP. Call it before the first OpenMP region.
inline int configure_threads(int requested) {
    int actual = 1;
#ifdef HEAT2D_USE_OPENMP
    omp_set_dynamic(0);
    if (requested > 0) {
        omp_set_num_threads(requested);
    } else if (!std::getenv("OMP_NUM_THREADS")) {
        omp_set_num_threads(1);
    }
#pragma omp parallel
    {
#pragma omp single
        actual = omp_get_num_threads();
    }
#else
    if (requested > 1) throw std::invalid_argument("this executable was built without OpenMP");
#endif
    return actual;
}

inline double sum_over_ranks(double value, MPI_Comm comm) {
    double total = 0.0;
    mpi_check(MPI_Allreduce(&value, &total, 1, MPI_DOUBLE, MPI_SUM, comm), "global sum");
    return total;
}

inline double max_over_ranks(double value, MPI_Comm comm) {
    double result = 0.0;
    mpi_check(MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, comm), "global maximum");
    return result;
}

inline double min_over_ranks(double value, MPI_Comm comm) {
    double result = 0.0;
    mpi_check(MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MIN, comm), "global minimum");
    return result;
}

inline int min_over_ranks(int value, MPI_Comm comm) {
    int result = 0;
    mpi_check(MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, comm), "global minimum");
    return result;
}

// File output happens on rank 0 only. If it fails there, the other ranks
// would wait forever in the next collective call, and unwinding the stack
// first would call MPI_Comm_free (a collective call) while they wait. So stop
// the whole job on the spot, with the message main() would have printed.
template <typename Action>
void root_io(const std::string& program, const ProcessLayout& layout, Action action) {
    if (layout.rank != 0) return;
    try {
        action();
    } catch (const std::exception& e) {
        std::cerr << program << " rank " << layout.rank << ": " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

// Every rank must run the same number of threads, or the timings and the
// thread count in the summary would mean nothing.
inline void check_same_threads(int threads, MPI_Comm comm) {
    int lowest = 0;
    int highest = 0;
    mpi_check(MPI_Allreduce(&threads, &lowest, 1, MPI_INT, MPI_MIN, comm), "thread count minimum");
    mpi_check(MPI_Allreduce(&threads, &highest, 1, MPI_INT, MPI_MAX, comm), "thread count maximum");
    if (lowest != highest) throw std::runtime_error("all ranks must use the same number of threads");
}

// The timings count one thread per physical core. Warn when a rank has fewer
// cores than threads, or when the threads are not bound one to each core.
inline void warn_about_binding(const std::string& program, int cores, Binding binding, int threads) {
    std::string problem;
    if (cores > 0 && cores < threads) {
        problem = "a rank may run on only " + std::to_string(cores) + (cores == 1 ? " core" : " cores") + " for its " +
                  std::to_string(threads) + " threads, so threads share cores";
    } else if (binding == Binding::shared) {
        problem = "two threads are bound to the same physical core (for example its two hardware threads)";
    } else if (binding == Binding::unpinned) {
        problem = "some threads are not bound to a single core, so the operating system decides where they run";
    } else {
        return;
    }
    std::cerr << program << ": warning: " << problem << ". Bind one thread to each physical core "
              << "(hpc/README.md, \"Binding\").\n";
}

// The first line of the MPI library's version string.
inline std::string mpi_library() {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int length = 0;
    mpi_check(MPI_Get_library_version(version, &length), "MPI_Get_library_version");
    std::string text(version, static_cast<std::size_t>(length));
    const std::size_t end = text.find_first_of("\n,");
    return end == std::string::npos ? text : text.substr(0, end);
}

// The strict interior for overlap (src/README.md, "Strict interior under
// overlap"): cells at least `halo` away from every edge shared with another
// rank. Cells next to a physical wall can be interior, because their ghosts
// are filled locally before the update. The 3D pigment model splits only x
// and y, so the same rows and columns are the strict interior of every layer.
struct Interior {
    int i_lo;
    int i_hi;
    int j_lo;
    int j_hi;
    bool empty() const { return i_lo > i_hi || j_lo > j_hi; }
};

inline Interior strict_interior(const ProcessLayout& layout) {
    const Block& b = layout.block;
    return Interior{layout.north != MPI_PROC_NULL ? 1 + halo : 1, layout.south != MPI_PROC_NULL ? b.ny - halo : b.ny,
                    layout.west != MPI_PROC_NULL ? 1 + halo : 1, layout.east != MPI_PROC_NULL ? b.nx - halo : b.nx};
}
