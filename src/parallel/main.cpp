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
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int exit_code = 0;
    if (wants_help(argc, argv)) {
        if (rank == 0) std::cout << usage("heat_parallel");
    } else {
        // Every rank parses the same arguments, so invalid input fails on all
        // ranks alike and they all reach MPI_Finalize.
        try {
            if (provided < MPI_THREAD_FUNNELED) {
                throw std::runtime_error("the MPI library does not provide MPI_THREAD_FUNNELED");
            }
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
        } catch (const std::exception& e) {
            if (rank == 0) {
                std::cerr << "heat_parallel: error: " << e.what() << "\n"
                          << "Run heat_parallel --help for the options.\n";
            }
            exit_code = 2;
        }
    }
    MPI_Finalize();
    return exit_code;
}
