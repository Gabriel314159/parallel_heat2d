#include "common/config.hpp"
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace heat2d {
namespace {
int integer(const std::string& value) {
    std::size_t used = 0;
    const long long result = std::stoll(value, &used);
    if (used != value.size() || result < 0 || result > std::numeric_limits<int>::max()-2)
        throw std::invalid_argument("Invalid nonnegative integer: " + value);
    return static_cast<int>(result);
}
double number(const std::string& value) {
    std::size_t used = 0;
    const double result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid finite number: " + value);
    return result;
}
}
Config parse_config(int argc, char** argv) {
    Config c;
    std::set<std::string> seen;
    for (int k=1; k<argc; ++k) {
        const std::string key = argv[k];
        if (key == "--help" || key == "-h") { c.help=true; return c; }
        if (!seen.insert(key).second) throw std::invalid_argument("Duplicate option: " + key);
        if (k+1 >= argc) throw std::invalid_argument("Missing value for " + key);
        const std::string value = argv[++k];
        if (key == "--nx") c.nx=integer(value);
        else if (key == "--ny") c.ny=integer(value);
        else if (key == "--steps") c.steps=integer(value);
        else if (key == "--lx") c.lx=number(value);
        else if (key == "--ly") c.ly=number(value);
        else if (key == "--alpha") c.alpha=number(value);
        else if (key == "--cfl") c.cfl=number(value);
        else if (key == "--dt") c.dt=number(value);
        else if (key == "--initial-temperature") c.initial=number(value);
        else if (key == "--top") c.top=number(value);
        else if (key == "--bottom") c.bottom=number(value);
        else if (key == "--left") c.left=number(value);
        else if (key == "--right") c.right=number(value);
        else if (key == "--amplitude") c.amplitude=number(value);
        else if (key == "--case") c.problem=value;
        else if (key == "--decomp") c.decomposition=value;
        else if (key == "--exchange") c.exchange=value;
        else if (key == "--py") c.py=integer(value);
        else if (key == "--px") c.px=integer(value);
        else if (key == "--threads") c.threads=integer(value);
        else if (key == "--output") c.output=value;
        else throw std::invalid_argument("Unknown option: " + key);
    }
    if (c.nx < 1 || c.ny < 1 || c.lx <= 0 || c.ly <= 0 || c.alpha <= 0)
        throw std::invalid_argument("Grid sizes, plate dimensions and alpha must be positive.");
    if (!(c.cfl > 0 && c.cfl <= 1)) throw std::invalid_argument("cfl must be in (0,1].");
    if (seen.count("--dt") && c.dt <= 0) throw std::invalid_argument("dt must be positive.");
    if (seen.count("--threads") && c.threads < 1) throw std::invalid_argument("threads must be positive.");
    if (c.decomposition != "rows" && c.decomposition != "2d")
        throw std::invalid_argument("decomp must be rows or 2d.");
    if (c.exchange != "blocking" && c.exchange != "overlap")
        throw std::invalid_argument("exchange must be blocking or overlap.");
    const std::string edges[] = {"--top", "--bottom", "--left", "--right"};
    double* values[] = {&c.top, &c.bottom, &c.left, &c.right};
    if (c.problem == "constant" || c.problem == "sine_decay") {
        const double expected = c.problem == "constant" ? c.initial : 0.0;
        for (int i=0; i<4; ++i) {
            if (!seen.count(edges[i])) *values[i]=expected;
            if (*values[i] != expected)
                throw std::invalid_argument(c.problem + " requires consistent fixed boundaries.");
        }
    } else if (c.problem != "heated_top") {
        throw std::invalid_argument("case must be heated_top, constant or sine_decay.");
    }
    c.dx=c.lx/(c.nx+1.0); c.dy=c.ly/(c.ny+1.0);
    if (c.dt == 0) c.dt=c.cfl/(2*c.alpha*(1/(c.dx*c.dx)+1/(c.dy*c.dy)));
    c.rx=c.alpha*c.dt/(c.dx*c.dx); c.ry=c.alpha*c.dt/(c.dy*c.dy);
    if (!(c.dt>0) || !std::isfinite(c.rx) || !std::isfinite(c.ry) ||
        !std::isfinite(c.steps*c.dt) || c.rx+c.ry > 0.5+32*std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("Unstable or unrepresentable time step: require rx+ry <= 0.5.");
    return c;
}
std::string usage() {
    return "Heat2D: explicit 2D heat diffusion\n"
           "  --nx N --ny N --steps N       Interior grid and time steps (steps may be zero)\n"
           "  --lx L --ly L --alpha A      Geometry and positive diffusivity\n"
           "  --cfl F | --dt DT            Automatic safety fraction or explicit time step\n"
           "  --case heated_top|constant|sine_decay\n"
           "  --initial-temperature T --top T --bottom T --left T --right T\n"
           "  --amplitude A                Sine initial-condition amplitude\n"
           "  --output PATH                Final CSV, including physical boundaries\n"
           "Parallel options:\n"
           "  --decomp rows|2d --exchange blocking|overlap\n"
           "  --py N --px N                Optional Cartesian process-grid dimensions\n"
           "  --threads N                  Override OMP_NUM_THREADS (default: env, then 1)\n"
           "MPI ranks are set by the MPI launcher. Results are printed as one JSON line.\n";
}
}
