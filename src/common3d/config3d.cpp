// Function bodies for common3d/config3d.hpp.

#include "common3d/config3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Replace the first occurrence of `from` in `text`, if any.
void replace_once(std::string& text, const std::string& from, const std::string& to) {
    const std::string::size_type at = text.find(from);
    if (at != std::string::npos) text.replace(at, from.size(), to);
}

// limit / rate, or infinity when the rate is zero.
double allowed(double limit, double rate) {
    return rate > 0 ? limit / rate : std::numeric_limits<double>::infinity();
}

struct Rates3D {
    double advection;  // a
    double diffusion;  // d, with the vertical term
    double viscosity;  // n, horizontal only
};

Rates3D rates_3d(const Config3D& c, double u_max, double v_max) {
    const double dx = c.base.lx / c.base.nx;
    const double dy = c.base.ly / c.base.ny;
    const double dz = dz_of(c);
    const double horizontal = 1 / (dx * dx) + 1 / (dy * dy);
    return Rates3D{u_max / dx + v_max / dy, 2 * c.base.kappa * (horizontal + 1 / (dz * dz)),
                   c.base.viscosity * horizontal};
}

}  // namespace

std::string usage_3d(const std::string& program) {
    std::string text = usage(program);
    replace_once(text, "Usage: " + program + " [options]\n",
                 "Usage: " + program + " [options]\n"
                 "\n"
                 "3D pigment in the 2D wind-driven flow (docs/design_3d.md, \"Model A\"). The pigment fills\n"
                 "a basin lx x ly x lz on nx x ny x nz cells. The flow is the 2D flow of heat_serial, the same\n"
                 "in every layer, with no vertical velocity; the pigment also diffuses vertically.\n"
                 "\n"
                 "Depth\n"
                 "  --nz N                  layers (default 32); layer 1 is at the surface\n"
                 "  --lz L                  depth (default 0.125)\n"
                 "  --blob-z Z|column       depth of the pigment_drop centre (default lz / 2); 'column' puts\n"
                 "                          the 2D drop in every layer, so each layer evolves as the 2D model\n"
                 "  --freeze-flow           advance the flow during the spin-up only; afterwards the pigment\n"
                 "                          moves in that fixed flow, with no Poisson solves\n"
                 "  --output PATH           the whole 3D field: CSV of nz x ny rows of nx values (layer by\n"
                 "                          layer, surface first) or .npy of shape (nz, ny, nx)\n"
                 "  Snapshots hold the surface layer (field_*.npy, with psi and the vorticity) and a vertical\n"
                 "  section through the row of the blob centre (section_*.npy, nz x nx).\n"
                 "\n"
                 "The other options are those of the 2D solver:\n");
    replace_once(text, "Parallel execution (heat_parallel only)", "Parallel execution (heat3d_parallel only)");
    replace_once(text, "  --decomp rows|2d        split rows only, or a 2D process grid (default rows)\n",
                 "  --decomp slabs|pencils  split y only, or x and y (default slabs); every rank keeps whole\n"
                 "                          columns of water, so z is never split\n");
    replace_once(text, "process grid rows and columns for --decomp 2d", "process grid rows and columns for --decomp pencils");
    return text;
}

std::string decomp_name_3d(Decomp d) {
    return d == Decomp::rows ? "slabs" : "pencils";
}

Config3D parse_args_3d(int argc, char** argv) {
    Config3D c;
    std::vector<std::string> rest;  // the 2D options, in their order
    rest.emplace_back(argc > 0 ? argv[0] : "heat3d");
    auto value_after = [&](int& k) -> std::string {
        const std::string flag = argv[k];
        if (k + 1 >= argc) throw std::invalid_argument("missing value after " + flag);
        ++k;
        return argv[k];
    };
    for (int k = 1; k < argc; ++k) {
        const std::string flag = argv[k];
        if (flag == "--nz") {
            c.nz = std::stoi(value_after(k));
        } else if (flag == "--lz") {
            c.lz = std::stod(value_after(k));
        } else if (flag == "--blob-z") {
            const std::string value = value_after(k);
            if (value == "column") {
                c.blob_column = true;
                c.blob_z.reset();
            } else {
                c.blob_column = false;
                c.blob_z = std::stod(value);
            }
        } else if (flag == "--freeze-flow") {
            c.freeze_flow = true;
        } else if (flag == "--decomp") {
            const std::string value = value_after(k);
            rest.emplace_back("--decomp");
            if (value == "slabs") {
                rest.emplace_back("rows");
            } else if (value == "pencils") {
                rest.emplace_back("2d");
            } else if (value == "blocks") {
                throw std::invalid_argument(
                    "--decomp blocks would split the depth, which needs the 3D flow (docs/design_3d.md, model B); "
                    "use slabs or pencils");
            } else {
                throw std::invalid_argument("unknown decomposition: " + value + " (use slabs or pencils)");
            }
        } else {
            rest.push_back(flag);
        }
    }
    std::vector<char*> words;
    for (std::string& word : rest) words.push_back(word.data());
    c.base = parse_args(static_cast<int>(words.size()), words.data());
    return c;
}

void validate_3d(const Config3D& c) {
    validate(c.base);
    if (c.nz < 1) throw std::invalid_argument("nz must be at least 1");
    if (!(c.lz > 0) || !std::isfinite(c.lz)) throw std::invalid_argument("lz must be positive");
    if (c.blob_z && !std::isfinite(*c.blob_z)) throw std::invalid_argument("blob_z must be a number or 'column'");
}

double blob_depth(const Config3D& c) {
    return c.blob_z.value_or(0.5 * c.lz);
}

double stable_dt_3d(const Config3D& c, double u_max, double v_max) {
    const Rates3D r = rates_3d(c, u_max, v_max);
    return std::min({allowed(1.0, r.advection + r.diffusion), allowed(0.5, r.advection), allowed(0.1, r.viscosity)});
}

void check_time_step_3d(const Config3D& c, double dt, double u_max, double v_max) {
    const Rates3D r = rates_3d(c, u_max, v_max);
    std::ostringstream problem;
    if (dt * (r.advection + r.diffusion) > 1.0) {
        problem << "pigment limit dt (a + d) = " << dt * (r.advection + r.diffusion) << " > 1";
    } else if (dt * r.advection > 0.5) {
        problem << "vorticity advection limit dt a = " << dt * r.advection << " > 0.5";
    } else if (dt * r.viscosity > 0.1) {
        problem << "viscosity limit dt nu (1/dx^2 + 1/dy^2) = " << dt * r.viscosity << " > 0.1";
    } else {
        return;
    }
    problem << " at measured speeds u_max = " << u_max << ", v_max = " << v_max
            << ". Rerun with a larger --max-speed or a smaller --dt.";
    throw std::runtime_error(problem.str());
}

double time_step_3d(const Config3D& c) {
    const double largest = stable_dt_3d(c, c.base.max_speed, c.base.max_speed);
    if (c.base.dt) {
        if (*c.base.dt > largest) {
            throw std::invalid_argument("dt = " + std::to_string(*c.base.dt) + " is unstable for --max-speed " +
                                        std::to_string(c.base.max_speed) + " and dz = " + std::to_string(dz_of(c)) +
                                        "; the limit is " + std::to_string(largest));
        }
        return *c.base.dt;
    }
    return c.base.cfl * largest;
}
