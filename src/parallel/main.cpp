// heat_parallel: the MPI + OpenMP solver. Run  heat_parallel --help  for the options.
//
//   OMP_NUM_THREADS=2 mpiexec -n 4 ./build/heat_parallel --decomp 2d --exchange overlap

#include <mpi.h>

#include <exception>
#include <iostream>
#include <stdexcept>

#include "common/config.hpp"
#include "common/io.hpp"
#include "parallel/solver.hpp"

int main(int argc, char** argv) {
    // FUNNELED: threads exist, but only the thread that initialized MPI calls it.
    int provided = 0;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS) {
        std::cerr << "heat_parallel: MPI initialization failed.\n";
        return 1;
    }
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // A failing MPI call now returns an error code, which mpi_check turns into
    // an exception, instead of aborting on the spot.
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

    try {
        if (provided < MPI_THREAD_FUNNELED) {
            throw std::runtime_error("the MPI library does not provide MPI_THREAD_FUNNELED");
        }
        if (wants_help(argc, argv)) {
            if (world_rank == 0) std::cout << usage("heat_parallel");
        } else {
            const Config cfg = parse_args(argc, argv);
            validate(cfg);

            const ParallelResult result = run_parallel(cfg, MPI_COMM_WORLD);
            if (result.is_root) {
                print_summary(std::cout, cfg, result.summary);
                if (!cfg.output.empty()) {
                    write_field(cfg.output, result.field, cfg.ny, cfg.nx);
                    write_text_file(replace_extension(cfg.output, ".json"), summary_json(cfg, result.summary) + "\n");
                }
            }
        }
    } catch (const std::exception& e) {
        // A failure on any rank stops all ranks. Otherwise the others would
        // wait forever in a collective call for the rank that has failed.
        std::cerr << "heat_parallel rank " << world_rank << ": " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
    return 0;
}
