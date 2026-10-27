#pragma once
#include "common/config.hpp"
#include "common/grid.hpp"

namespace heat2d {
double initial_value(const Config& c, int global_row, int global_column);
void initialize(Grid& g, const Config& c, int offset_y=0, int offset_x=0);
double analytical_value(const Config& c, int global_row, int global_column, double time);
double analytical_error(const Grid& g, const Config& c, int offset_y=0, int offset_x=0);
}
