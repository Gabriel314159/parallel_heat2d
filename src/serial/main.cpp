#include "common/config.hpp"
#include "common/io.hpp"
#include "common/physics.hpp"
#include "serial/solver.hpp"
#include <exception>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        const auto c=heat2d::parse_config(argc,argv);
        if (c.help) { std::cout << heat2d::usage(); return 0; }
        if (c.threads>1) throw std::invalid_argument("heat_serial always uses one thread.");
        heat2d::Grid current(c.nx,c.ny), next(c.nx,c.ny);
        heat2d::initialize(current,c); heat2d::initialize(next,c);
        heat2d::RunInfo result;
        result.elapsed=heat2d::solve_serial(current,next,c);
        if (c.problem=="sine_decay") result.analytical_max_abs_error=heat2d::analytical_error(current,c);
        if (!c.output.empty()) heat2d::write_field(c.output,current);
        heat2d::print_result(c,result);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "heat_serial: " << e.what() << '\n'; return 1;
    }
}
