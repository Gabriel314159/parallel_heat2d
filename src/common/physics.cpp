#include "common/physics.hpp"
#include <algorithm>
#include <cmath>

namespace heat2d {
double initial_value(const Config& c, int i, int j) {
    // Left/right values take precedence at corners in every implementation.
    if (j==0) return c.left;
    if (j==c.nx+1) return c.right;
    if (i==0) return c.top;
    if (i==c.ny+1) return c.bottom;
    if (c.problem=="sine_decay") {
        const double pi=std::acos(-1.0);
        return c.amplitude*std::sin(pi*j/(c.nx+1.0))*std::sin(pi*i/(c.ny+1.0));
    }
    return c.initial;
}
void initialize(Grid& g, const Config& c, int oy, int ox) {
    for (int i=0; i<=g.ny+1; ++i)
        for (int j=0; j<=g.nx+1; ++j) g(i,j)=initial_value(c,oy+i,ox+j);
}
double analytical_value(const Config& c, int i, int j, double time) {
    const double pi=std::acos(-1.0);
    return initial_value(c,i,j)*std::exp(-c.alpha*pi*pi*(1/(c.lx*c.lx)+1/(c.ly*c.ly))*time);
}
double analytical_error(const Grid& g, const Config& c, int oy, int ox) {
    double error=0;
    for (int i=1; i<=g.ny; ++i)
        for (int j=1; j<=g.nx; ++j)
            error=std::max(error,std::abs(g(i,j)-analytical_value(c,oy+i,ox+j,c.steps*c.dt)));
    return error;
}
}
