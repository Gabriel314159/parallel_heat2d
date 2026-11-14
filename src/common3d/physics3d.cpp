// Function bodies for common3d/physics3d.hpp.

#include "common3d/physics3d.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace {

constexpr double pi = M_PI;

}  // namespace

double initial_value_3d(const Config3D& c, double x, double y, double z) {
    const Config& b = c.base;
    switch (b.initial_case) {
        case InitialCase::pigment_drop: {
            if (c.blob_column) return initial_value(b, x, y);
            const double rx = x - b.blob_x;
            const double ry = y - b.blob_y;
            const double rz = z - blob_depth(c);
            const double r2 = rx * rx + ry * ry + rz * rz;
            return b.background + b.amplitude * std::exp(-r2 / (2 * b.blob_width * b.blob_width));
        }
        case InitialCase::constant:
            return b.background;
        case InitialCase::cosine_decay:
            return b.amplitude * std::cos(pi * x / b.lx) * std::cos(pi * y / b.ly) * std::cos(pi * z / c.lz);
    }
    return 0.0;
}

void set_initial_condition_3d(Field3D& field, const Config3D& c, const Block& b) {
    assert(field.nz() == c.nz && field.ny() == b.ny && field.nx() == b.nx);
    const double dx = c.base.lx / c.base.nx;
    const double dy = c.base.ly / c.base.ny;
    const double dz = dz_of(c);
    for (int k = 1; k <= c.nz; ++k) {
        const double z = cell_centre_z(k, dz, c.lz);
        for (int i = 1; i <= b.ny; ++i) {
            const double y = cell_centre_y(b.row_offset + i, dy, c.base.ly);
            for (int j = 1; j <= b.nx; ++j) {
                const double x = cell_centre_x(b.col_offset + j, dx);
                field(k, i, j) = initial_value_3d(c, x, y, z);
            }
        }
    }
}

void fill_wall_ghosts_3d(Field3D& field, const Block& b) {
    static_assert(halo == 1, "the wall ghosts below mirror one layer");
    const int nz = field.nz();
    const int ny = field.ny();
    const int nx = field.nx();
    // Side walls, layer by layer: ghost = the interior cell across the wall.
    for (int k = 1; k <= nz; ++k) {
        if (b.wall_top) {
            for (int j = 1; j <= nx; ++j) field(k, 0, j) = field(k, 1, j);
        }
        if (b.wall_bottom) {
            for (int j = 1; j <= nx; ++j) field(k, ny + 1, j) = field(k, ny, j);
        }
        if (b.wall_left) {
            for (int i = 1; i <= ny; ++i) field(k, i, 0) = field(k, i, 1);
        }
        if (b.wall_right) {
            for (int i = 1; i <= ny; ++i) field(k, i, nx + 1) = field(k, i, nx);
        }
    }
    // The surface (above layer 1) and the bottom (below layer nz).
    for (int i = 1; i <= ny; ++i) {
        const double* top = field.row(1, i);
        const double* bottom = field.row(nz, i);
        double* above = field.row(0, i);
        double* below = field.row(nz + 1, i);
        for (int j = 1; j <= nx; ++j) {
            above[j] = top[j];
            below[j] = bottom[j];
        }
    }
}

double local_mass_3d(const Field3D& field, double dx, double dy, double dz) {
    double sum = 0.0;
    for (int k = 1; k <= field.nz(); ++k) {
        for (int i = 1; i <= field.ny(); ++i) {
            const double* row = field.row(k, i);
            for (int j = 1; j <= field.nx(); ++j) sum += row[j];
        }
    }
    return sum * dx * dy * dz;
}

double local_abs_mass_3d(const Field3D& field, double dx, double dy, double dz) {
    double sum = 0.0;
    for (int k = 1; k <= field.nz(); ++k) {
        for (int i = 1; i <= field.ny(); ++i) {
            const double* row = field.row(k, i);
            for (int j = 1; j <= field.nx(); ++j) sum += std::abs(row[j]);
        }
    }
    return sum * dx * dy * dz;
}

double cosine_decay_continuous_3d(const Config3D& c, double x, double y, double z, double t) {
    const Config& b = c.base;
    const double rate = b.kappa * pi * pi * (1 / (b.lx * b.lx) + 1 / (b.ly * b.ly) + 1 / (c.lz * c.lz));
    return b.amplitude * std::cos(pi * x / b.lx) * std::cos(pi * y / b.ly) * std::cos(pi * z / c.lz) *
           std::exp(-rate * t);
}

double cosine_decay_discrete_factor_3d(const Config3D& c, double dt, int steps) {
    const Config& b = c.base;
    const double dx = b.lx / b.nx;
    const double dy = b.ly / b.ny;
    const double dz = dz_of(c);
    const double sx = std::sin(pi * dx / (2 * b.lx));
    const double sy = std::sin(pi * dy / (2 * b.ly));
    const double sz = std::sin(pi * dz / (2 * c.lz));
    const double lambda = 4 * b.kappa / (dx * dx) * sx * sx + 4 * b.kappa / (dy * dy) * sy * sy +
                          4 * b.kappa / (dz * dz) * sz * sz;
    return std::pow(1 - lambda * dt, steps);
}

CosineErrors cosine_decay_errors_3d(const Field3D& field, const Config3D& c, const Block& b, double dt, int steps) {
    const double dx = c.base.lx / c.base.nx;
    const double dy = c.base.ly / c.base.ny;
    const double dz = dz_of(c);
    const double factor = cosine_decay_discrete_factor_3d(c, dt, steps);
    const double t = steps * dt;
    const double scale = c.base.amplitude != 0.0 ? std::abs(c.base.amplitude) : 1.0;

    CosineErrors e;
    for (int k = 1; k <= c.nz; ++k) {
        const double z = cell_centre_z(k, dz, c.lz);
        for (int i = 1; i <= b.ny; ++i) {
            const double y = cell_centre_y(b.row_offset + i, dy, c.base.ly);
            for (int j = 1; j <= b.nx; ++j) {
                const double x = cell_centre_x(b.col_offset + j, dx);
                const double value = field(k, i, j);
                e.discrete = std::max(e.discrete, std::abs(value - factor * initial_value_3d(c, x, y, z)));
                e.continuous = std::max(e.continuous, std::abs(value - cosine_decay_continuous_3d(c, x, y, z, t)));
            }
        }
    }
    e.discrete /= scale;
    e.continuous /= scale;
    return e;
}
