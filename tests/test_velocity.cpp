// Unit test for src/common/velocity.cpp. Built and run by CTest: ctest --test-dir build

#include "common/velocity.hpp"

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

    // FlowField (the fast fill used by the solvers) must agree exactly with the
    // per-face functions above on interior faces, and be exactly 0 on walls.
    // Check the whole grid and a sub-block with offsets, as a parallel rank sees it.
    {
        const double t = 2.5;
        Block whole;
        whole.ny = ny;
        whole.nx = nx;
        Block part;  // rows 11..20, columns 6..17 of the global grid, walls nowhere
        part.ny = 10;
        part.nx = 12;
        part.row_offset = 10;
        part.col_offset = 5;
        part.wall_top = part.wall_bottom = part.wall_left = part.wall_right = false;

        for (const Block& b : {whole, part}) {
            FlowField flow(p, b, dx, dy);
            flow.set_time(t);
            FaceVelocities vel(b.ny, b.nx);
            flow.fill_all(vel);

            bool exact = true;
            for (int i = 1; i <= b.ny; ++i) {
                for (int j = 0; j <= b.nx; ++j) {
                    const bool wall = (j == 0 && b.wall_left) || (j == b.nx && b.wall_right);
                    const double expected =
                        wall ? 0.0 : face_velocity_east(p, b.row_offset + i, b.col_offset + j, dx, dy, t);
                    exact = exact && vel.east(i, j) == expected;
                }
            }
            for (int i = 1; i <= b.ny + 1; ++i) {
                for (int j = 1; j <= b.nx; ++j) {
                    const bool wall = (i == 1 && b.wall_top) || (i == b.ny + 1 && b.wall_bottom);
                    const double expected =
                        wall ? 0.0 : face_velocity_north(p, b.row_offset + i, b.col_offset + j, dx, dy, t);
                    exact = exact && vel.north(i, j) == expected;
                }
            }
            std::cout << "FlowField on " << b.ny << " x " << b.nx << " block at (" << b.row_offset << ", "
                      << b.col_offset << "): " << (exact ? "matches exactly" : "MISMATCH") << "\n";
            ok = ok && exact;
        }
        FlowField steady_flow(FlowParams{0.1, 0.0, p.omega, p.lx, p.ly}, whole, dx, dy);
        FlowField no_flow(FlowParams{0.0, 0.25, p.omega, p.lx, p.ly}, whole, dx, dy);
        ok = ok && steady_flow.steady() && no_flow.steady() && !FlowField(p, whole, dx, dy).steady();
    }

    std::cout << "\n" << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 1;
}
