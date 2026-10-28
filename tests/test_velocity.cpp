// Unit test for src/common/velocity.cpp. Built and run by CTest: ctest --test-dir build

#include "common/velocity.hpp"

#include "common/corners.hpp"
#include "common/flow.hpp"
#include "common/grid.hpp"

#include <algorithm>  // std::max
#include <cmath>
#include <iostream>

int main() {
    // Defaults from config/defaults.json, on a deliberately uneven small grid.
    FlowParams p{0.1, 0.25, 2.0 * M_PI / 10.0, 1.0, 1.0};
    const int nx = 37;
    const int ny = 29;
    const double dx = p.lx / nx;
    const double dy = p.ly / ny;

    const double u_bound = max_speed_u(p);
    const double v_bound = max_speed_v(p);
    std::cout << "u_max bound = " << u_bound << "   (expect 0.314159)\n";
    std::cout << "v_max bound = " << v_bound << "   (expect 0.942478)\n";
    std::cout << "psi(0.3, 0.4, 1.7) = " << streamfunction(p, 0.3, 0.4, 1.7)
              << "   (expect 0.0918113)\n\n";

    bool ok = u_bound > 0.0 && v_bound > 0.0;

    // Check several times within one flow period.
    for (double t : {0.0, 1.3, 2.5, 7.9}) {
        double max_div = 0.0;
        double max_u = 0.0;
        double max_v = 0.0;
        double max_wall = 0.0;

        for (int i = 1; i <= ny; ++i) {
            for (int j = 1; j <= nx; ++j) {
                const double ue = face_velocity_east(p, i, j, dx, dy, t);
                const double uw = face_velocity_east(p, i, j - 1, dx, dy, t);
                const double vn = face_velocity_north(p, i, j, dx, dy, t);
                const double vs = face_velocity_north(p, i + 1, j, dx, dy, t);

                // Discrete divergence of cell (i, j): should be rounding-level.
                const double div = (ue - uw) / dx + (vn - vs) / dy;
                max_div = std::max(max_div, std::abs(div));
                max_u = std::max(max_u, std::abs(ue));
                max_v = std::max(max_v, std::abs(vn));

                // Faces lying on the walls: normal velocity should vanish.
                if (j == 1) max_wall = std::max(max_wall, std::abs(uw));
                if (j == nx) max_wall = std::max(max_wall, std::abs(ue));
                if (i == 1) max_wall = std::max(max_wall, std::abs(vn));
                if (i == ny) max_wall = std::max(max_wall, std::abs(vs));
            }
        }

        std::cout << "t = " << t << "\n"
                  << "  max |divergence|      = " << max_div << "\n"
                  << "  max |U_e| / u bound   = " << max_u / u_bound << "\n"
                  << "  max |V_n| / v bound   = " << max_v / v_bound << "\n"
                  << "  max |wall velocity|   = " << max_wall << "\n";

        ok = ok && max_div < 1e-12 && max_u <= u_bound && max_v <= v_bound
             && max_wall < 1e-12 && max_v > 0.0;
    }

    // faces_from_psi (used by the solvers) applies the same face formulas to psi
    // stored at the corners. With psi sampled from the analytic streamfunction,
    // it must reproduce the per-face functions above exactly, on the whole grid
    // and on a sub-block with offsets, as a parallel rank sees it.
    {
        const double t = 2.5;
        Block whole;
        whole.ny = ny;
        whole.nx = nx;
        Block part;  // rows 11..20, columns 6..17 of the global grid
        part.ny = 10;
        part.nx = 12;
        part.row_offset = 10;
        part.col_offset = 5;
        part.wall_top = part.wall_bottom = part.wall_left = part.wall_right = false;

        for (const Block& b : {whole, part}) {
            Field psi(b.ny, b.nx);
            for (int i = 0; i <= b.ny; ++i) {
                for (int j = 0; j <= b.nx; ++j) {
                    psi(i, j) = streamfunction(p, corner_x(b.col_offset + j, dx), corner_y(b.row_offset + i, dy, p.ly), t);
                }
            }
            FaceVelocities vel(b.ny, b.nx);
            faces_from_psi(psi, vel, dx, dy);

            bool exact = true;
            for (int i = 1; i <= b.ny; ++i) {
                for (int j = 0; j <= b.nx; ++j) {
                    exact = exact && vel.east(i, j) == face_velocity_east(p, b.row_offset + i, b.col_offset + j, dx, dy, t);
                }
            }
            for (int i = 1; i <= b.ny + 1; ++i) {
                for (int j = 1; j <= b.nx; ++j) {
                    exact = exact && vel.north(i, j) == face_velocity_north(p, b.row_offset + i, b.col_offset + j, dx, dy, t);
                }
            }
            std::cout << "faces_from_psi on " << b.ny << " x " << b.nx << " block at (" << b.row_offset << ", "
                      << b.col_offset << "): " << (exact ? "matches exactly" : "MISMATCH") << "\n";
            ok = ok && exact;
        }
    }

    std::cout << "\n" << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 1;
}
