#pragma once
#include "common/config.hpp"
#include "common/grid.hpp"
namespace heat2d { double solve_serial(Grid& current, Grid& next, const Config& c); }
