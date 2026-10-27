#pragma once
#include "parallel/decomposition.hpp"
namespace heat2d {
double solve_parallel(Grid& current, Grid& next, const Config& c, const Domain& d);
int configure_threads(int requested);
}
