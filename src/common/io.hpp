#pragma once
#include "common/config.hpp"
#include "common/grid.hpp"
#include <string>

namespace heat2d {
struct RunInfo {
    std::string backend="serial", decomposition="none", exchange="none", mpi_library="none";
    int ranks=1, threads=1, nodes=1, py=1, px=1;
    double elapsed=0.0, analytical_max_abs_error=-1.0;
};
void write_field(const std::string& path, const Grid& g);
void print_result(const Config& c, const RunInfo& r);
}
