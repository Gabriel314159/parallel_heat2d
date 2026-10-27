#include "serial/solver.hpp"
#include "common/stencil.hpp"
#include <chrono>

namespace heat2d {
double solve_serial(Grid& current, Grid& next, const Config& c) {
    const auto start=std::chrono::steady_clock::now();
    for (int step=0; step<c.steps; ++step) {
        for (int i=1; i<=current.ny; ++i)
            for (int j=1; j<=current.nx; ++j)
                next(i,j)=update_point(current,i,j,c.rx,c.ry);
        current.swap_data(next);
    }
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
}
