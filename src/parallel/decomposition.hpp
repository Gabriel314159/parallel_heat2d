#pragma once
#include "common/config.hpp"
#include "common/grid.hpp"
#include <mpi.h>
#include <memory>
#include <stdexcept>
#include <string>

namespace heat2d {
inline void mpi_check(int code, const char* operation) {
    if (code != MPI_SUCCESS) {
        char text[MPI_MAX_ERROR_STRING]; int length=0;
        MPI_Error_string(code,text,&length);
        throw std::runtime_error(std::string(operation)+": "+std::string(text,length));
    }
}
class Domain {
public:
    MPI_Comm comm=MPI_COMM_NULL;
    int rank=0, size=1, py=1, px=1;
    int ny=0, nx=0, oy=0, ox=0;
    int north=MPI_PROC_NULL, south=MPI_PROC_NULL, west=MPI_PROC_NULL, east=MPI_PROC_NULL;
    Domain(const Config& c, MPI_Comm world);
    ~Domain();
    Domain(const Domain&)=delete;
    Domain& operator=(const Domain&)=delete;
};
std::unique_ptr<Grid> gather_global(const Grid& local, const Config& c, const Domain& d);
}
