// heat_serial: parse the command line, run the serial solver, report and
// write the requested output. Run  heat_serial --help  for the options.

#include <exception>
#include <iostream>
#include <stdexcept>

#include "common/config.hpp"
#include "common/io.hpp"
#include "serial/solver.hpp"

int main(int argc, char** argv) {
    if (wants_help(argc, argv)) {
        std::cout << usage("heat_serial");
        return 0;
    }
    try {
        const Config cfg = parse_args(argc, argv);
        validate(cfg);
        if (cfg.threads && *cfg.threads > 1) {
            throw std::invalid_argument("heat_serial always uses one thread");
        }
        if (has_parallel_options(cfg)) {
            std::cerr << "heat_serial: note: --decomp, --exchange, --py and --px are ignored (heat_parallel only)\n";
        }

        const SerialResult result = run_serial(cfg);
        print_summary(std::cout, cfg, result.summary);

        if (!cfg.output.empty()) {
            write_field(cfg.output, result.field, cfg.ny, cfg.nx);
            write_text_file(replace_extension(cfg.output, ".json"), summary_json(cfg, result.summary) + "\n");
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "heat_serial: error: " << e.what() << "\n"
                  << "Run heat_serial --help for the options.\n";
        return 2;
    }
}
