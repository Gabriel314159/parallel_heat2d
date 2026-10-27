#include "common/config.hpp"
#include "common/io.hpp"
#include "common/physics.hpp"
#include "parallel/solver.hpp"
#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    int provided=0;
    if (MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided)!=MPI_SUCCESS) {
        std::cerr << "MPI initialization failed.\n"; return 1;
    }
    int world_rank=0; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
    try {
        if (provided<MPI_THREAD_FUNNELED) throw std::runtime_error("MPI_THREAD_FUNNELED is required.");
        const auto c=heat2d::parse_config(argc,argv);
        if (c.help) {
            if (world_rank==0) std::cout << heat2d::usage();
        } else {
            heat2d::Domain d(c,MPI_COMM_WORLD);
            heat2d::Grid current(d.nx,d.ny), next(d.nx,d.ny);
            heat2d::initialize(current,c,d.oy,d.ox); heat2d::initialize(next,c,d.oy,d.ox);
            heat2d::RunInfo result;
            result.backend="parallel"; result.ranks=d.size;
            result.threads=heat2d::configure_threads(c.threads);
            int minimum=0, maximum=0;
            heat2d::mpi_check(MPI_Allreduce(&result.threads,&minimum,1,MPI_INT,MPI_MIN,d.comm),"threads min");
            heat2d::mpi_check(MPI_Allreduce(&result.threads,&maximum,1,MPI_INT,MPI_MAX,d.comm),"threads max");
            if (minimum!=maximum) throw std::runtime_error("All ranks must use the same number of threads.");
            result.py=d.py; result.px=d.px;
            result.decomposition=c.decomposition; result.exchange=c.exchange;
            MPI_Comm shared=MPI_COMM_NULL;
            heat2d::mpi_check(MPI_Comm_split_type(d.comm,MPI_COMM_TYPE_SHARED,d.rank,MPI_INFO_NULL,&shared),"shared communicator");
            int shared_rank=0; MPI_Comm_rank(shared,&shared_rank);
            int leader=shared_rank==0 ? 1:0;
            heat2d::mpi_check(MPI_Allreduce(&leader,&result.nodes,1,MPI_INT,MPI_SUM,d.comm),"node count");
            MPI_Comm_free(&shared);
            char library[MPI_MAX_LIBRARY_VERSION_STRING]; int length=0;
            MPI_Get_library_version(library,&length); result.mpi_library.assign(library,length);
            result.elapsed=heat2d::solve_parallel(current,next,c,d);
            if (c.problem=="sine_decay") {
                const double local=heat2d::analytical_error(current,c,d.oy,d.ox);
                heat2d::mpi_check(MPI_Allreduce(&local,&result.analytical_max_abs_error,1,MPI_DOUBLE,MPI_MAX,d.comm),"analytical error");
            }
            if (!c.output.empty()) {
                auto global=heat2d::gather_global(current,c,d);
                if (d.rank==0) heat2d::write_field(c.output,*global);
            }
            if (d.rank==0) heat2d::print_result(c,result);
        }
    } catch (const std::exception& e) {
        std::cerr << "heat_parallel rank " << world_rank << ": " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD,1); return 1;
    }
    MPI_Finalize(); return 0;
}
