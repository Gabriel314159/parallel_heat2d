// Corner fields: psi and zeta live on the cell corners.
// See docs/design.md, "Corner points".
//
// A corner field is an ordinary Field of the same shape as the pigment field.
// Entry (i, j) holds the south-east corner of cell (i, j), which is the corner
// with global indices (row_offset + i, col_offset + j). So the blocks, the
// decomposition and the halo exchange are reused unchanged.

#pragma once

#include "common/grid.hpp"

// Position of the corner with global indices (ig, jg).
inline double corner_x(int jg, double dx) { return jg * dx; }
inline double corner_y(int ig, double dy, double ly) { return ly - ig * dy; }

// The unknown corners of a block, as inclusive local index ranges. Corners on
// the bottom and right walls are stored in the block's last row and column
// (when the block lies on those walls) and always hold 0.
struct CornerRange {
    int i_first;
    int i_last;
    int j_first;
    int j_last;
    bool empty() const { return i_first > i_last || j_first > j_last; }
};

inline CornerRange unknown_corners(const Block& b) {
    return CornerRange{1, b.wall_bottom ? b.ny - 1 : b.ny, 1, b.wall_right ? b.nx - 1 : b.nx};
}

// Corners on the top and left walls are the ghost row and column of the blocks
// on those walls: set them, and the ghost corners next to them, to 0. Call this
// before a halo exchange with corners, so that zeros travel to the diagonal
// neighbours too.
inline void zero_wall_ghost_corners(Field& f, const Block& b) {
    for (int k = 1; k <= halo; ++k) {
        if (b.wall_top) {
            for (int j = 1 - halo; j <= f.nx() + halo; ++j) f(1 - k, j) = 0.0;
        }
        if (b.wall_left) {
            for (int i = 1 - halo; i <= f.ny() + halo; ++i) f(i, 1 - k) = 0.0;
        }
    }
}
