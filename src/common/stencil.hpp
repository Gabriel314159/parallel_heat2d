#pragma once
#include "common/grid.hpp"

namespace heat2d {
inline double update_point(const Grid& g, int i, int j, double rx, double ry) {
    const double center=g(i,j);
    return center + rx*(g(i,j-1)-2.0*center+g(i,j+1))
                  + ry*(g(i-1,j)-2.0*center+g(i+1,j));
}
}
