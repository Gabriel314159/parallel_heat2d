#include "common/io.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#ifndef _WIN32
#include <unistd.h>
#endif
#ifndef HEAT2D_BUILD_TYPE
#define HEAT2D_BUILD_TYPE "unknown"
#endif

namespace heat2d {
namespace {
std::string quoted(const std::string& s) {
    std::ostringstream out; out << '"';
    for (unsigned char ch:s) {
        if (ch=='"' || ch=='\\') out << '\\' << ch;
        else if (ch<32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(ch);
        else out << ch;
    }
    out << '"'; return out.str();
}
std::string hostname() {
#ifndef _WIN32
    char name[256]={};
    if (gethostname(name,sizeof(name)-1)==0) return name;
#endif
    const char* fallback=std::getenv("COMPUTERNAME");
    return fallback ? fallback : "unknown";
}
std::string compiler() {
#ifdef __VERSION__
    return __VERSION__;
#else
    return "unknown";
#endif
}
}
void write_field(const std::string& path, const Grid& g) {
    const std::filesystem::path file(path);
    if (!file.parent_path().empty()) std::filesystem::create_directories(file.parent_path());
    std::ofstream out(file);
    if (!out) throw std::runtime_error("Cannot open output: "+path);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int i=0; i<=g.ny+1; ++i) {
        for (int j=0; j<=g.nx+1; ++j) { if (j) out << ','; out << g(i,j); }
        out << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("Failed writing output: "+path);
}
void print_result(const Config& c, const RunInfo& r) {
    std::cout << std::setprecision(17)
      << "{\"backend\":" << quoted(r.backend)
      << ",\"nx\":" << c.nx << ",\"ny\":" << c.ny << ",\"steps\":" << c.steps
      << ",\"lx\":" << c.lx << ",\"ly\":" << c.ly << ",\"alpha\":" << c.alpha
      << ",\"dt\":" << c.dt << ",\"final_time\":" << c.steps*c.dt
      << ",\"rx\":" << c.rx << ",\"ry\":" << c.ry << ",\"case\":" << quoted(c.problem)
      << ",\"mpi_ranks\":" << r.ranks << ",\"omp_threads\":" << r.threads
      << ",\"total_cores\":" << r.ranks*r.threads << ",\"nodes\":" << r.nodes
      << ",\"py\":" << r.py << ",\"px\":" << r.px
      << ",\"decomposition\":" << quoted(r.decomposition) << ",\"exchange\":" << quoted(r.exchange)
      << ",\"elapsed_s\":" << r.elapsed << ",\"hostname\":" << quoted(hostname())
      << ",\"compiler\":" << quoted(compiler()) << ",\"build_type\":" << quoted(HEAT2D_BUILD_TYPE)
      << ",\"mpi_library\":" << quoted(r.mpi_library);
    if (r.analytical_max_abs_error>=0)
        std::cout << ",\"analytical_max_abs_error\":" << r.analytical_max_abs_error;
    std::cout << "}\n";
}
}
