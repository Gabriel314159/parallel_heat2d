// heat3d_parallel: the MPI + OpenMP 3D pigment solver (docs/design_3d.md,
// "Model A"). Run  heat3d_parallel --help  for the options.
//
//   OMP_NUM_THREADS=2 mpiexec -n 4 ./build/heat3d_parallel --nz 32 --decomp pencils --exchange overlap

#include <mpi.h>

#include <exception>
#include <iostream>
#include <stdexcept>

#include "common/config.hpp"
#include "common/io.hpp"
#include "common3d/config3d.hpp"
#include "common3d/io3d.hpp"
#include "parallel3d/solver.hpp"

int main(int argc, char** argv) {
    // FUNNELED: threads exist, but only the thread that initialized MPI calls it.
    int provided = 0;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS) {
        std::cerr << "heat3d_parallel: MPI initialization failed.\n";
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
            if (world_rank == 0) std::cout << usage_3d("heat3d_parallel");
        } else {
            const Config3D c3 = parse_args_3d(argc, argv);
            validate_3d(c3);
            const Config& cfg = c3.base;

            const Parallel3DResult result = run_parallel_3d(c3, MPI_COMM_WORLD);
            if (result.is_root) {
                print_summary(std::cout, cfg, result.summary);
                if (!cfg.output.empty()) {
                    write_field_3d(cfg.output, result.field, c3.nz, cfg.ny, cfg.nx);
                    write_text_file(replace_extension(cfg.output, ".json"), summary_json(cfg, result.summary) + "\n");
                }
            }
        }
    } catch (const std::exception& e) {
        // A failure on any rank stops all ranks. Otherwise the others would
        // wait forever in a collective call for the rank that has failed.
        std::cerr << "heat3d_parallel rank " << world_rank << ": " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
    return 0;
}
