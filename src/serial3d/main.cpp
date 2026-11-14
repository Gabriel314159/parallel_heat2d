// heat3d_serial: the 3D pigment in the 2D wind-driven flow, on one core
// (docs/design_3d.md, "Model A"). Run  heat3d_serial --help  for the options.
//
//   ./build/heat3d_serial --nx 64 --ny 64 --nz 16 --lz 0.25 --steps 200 --spinup 0 --output outputs/final3d.npy

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "common/config.hpp"
#include "common/io.hpp"
#include "common3d/config3d.hpp"
#include "common3d/io3d.hpp"
#include "serial3d/solver.hpp"

int main(int argc, char** argv) {
    const std::string name = serial3d_solver_name;
    if (wants_help(argc, argv)) {
        std::cout << usage_3d(name);
        return 0;
    }
    try {
        const Config3D c3 = parse_args_3d(argc, argv);
        validate_3d(c3);
        const Config& cfg = c3.base;
        if (cfg.threads && *cfg.threads > 1) {
            throw std::invalid_argument(name + " always uses one thread");
        }
        if (has_parallel_options(cfg)) {
            std::cerr << name << ": note: --decomp, --exchange, --py and --px are ignored (heat3d_parallel only)\n";
        }

        const Serial3DResult result = run_serial_3d(c3);
        print_summary(std::cout, cfg, result.summary);

        if (!cfg.output.empty()) {
            write_field_3d(cfg.output, result.field, c3.nz, cfg.ny, cfg.nx);
            write_text_file(replace_extension(cfg.output, ".json"), summary_json(cfg, result.summary) + "\n");
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << name << ": error: " << e.what() << "\n"
                  << "Run " << name << " --help for the options.\n";
        return 2;
    }
}
