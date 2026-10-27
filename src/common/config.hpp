#pragma once
#include <string>

namespace heat2d {
struct Config {
    int nx = 256, ny = 256, steps = 1000;
    double lx = 1.0, ly = 1.0, alpha = 1.0, cfl = 0.8, dt = 0.0;
    double initial = 20.0, top = 100.0, bottom = 20.0, left = 20.0, right = 20.0;
    double amplitude = 1.0, dx = 0.0, dy = 0.0, rx = 0.0, ry = 0.0;
    std::string problem = "heated_top", decomposition = "rows", exchange = "blocking";
    std::string output;
    int py = 0, px = 0, threads = 0;
    bool help = false;
};
Config parse_config(int argc, char** argv);
std::string usage();
}
