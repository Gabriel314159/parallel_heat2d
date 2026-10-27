// Function bodies for common/physics.hpp.

#include "common/physics.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

constexpr double pi = M_PI;

Block whole_domain(const Config& c) {
    Block b;
    b.ny = c.ny;
    b.nx = c.nx;
    return b;  // offsets 0, walls on all sides: the defaults of Block
}

double initial_value(const Config& c, double x, double y) {
    switch (c.initial_case) {
        case InitialCase::pigment_drop: {
            const double rx = x - c.blob_x;
            const double ry = y - c.blob_y;
            const double r2 = rx * rx + ry * ry;
            return c.background + c.amplitude * std::exp(-r2 / (2 * c.blob_width * c.blob_width));
        }
        case InitialCase::constant:
            return c.background;
        case InitialCase::cosine_decay:
            return c.amplitude * std::cos(pi * x / c.lx) * std::cos(pi * y / c.ly);
    }
    return 0.0;
}

void set_initial_condition(Field& field, const Config& c, const Block& b) {
    assert(field.ny() == b.ny && field.nx() == b.nx);
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;
    for (int i = 1; i <= b.ny; ++i) {
        const double y = cell_centre_y(b.row_offset + i, dy, c.ly);
        for (int j = 1; j <= b.nx; ++j) {
            const double x = cell_centre_x(b.col_offset + j, dx);
            field(i, j) = initial_value(c, x, y);
        }
    }
}

void fill_wall_ghosts(Field& field, const Block& b) {
    const int ny = field.ny();
    const int nx = field.nx();
    // Ghost layer k (1..halo) mirrors interior layer k across the wall face.
    for (int k = 1; k <= halo; ++k) {
        if (b.wall_top) {
            for (int j = 1; j <= nx; ++j) field(1 - k, j) = field(k, j);
        }
        if (b.wall_bottom) {
            for (int j = 1; j <= nx; ++j) field(ny + k, j) = field(ny + 1 - k, j);
        }
        if (b.wall_left) {
            for (int i = 1; i <= ny; ++i) field(i, 1 - k) = field(i, k);
        }
        if (b.wall_right) {
            for (int i = 1; i <= ny; ++i) field(i, nx + k) = field(i, nx + 1 - k);
        }
    }
}

double local_mass(const Field& field, double dx, double dy) {
    double sum = 0.0;
    for (int i = 1; i <= field.ny(); ++i) {
        const double* row = field.row(i);
        for (int j = 1; j <= field.nx(); ++j) sum += row[j];
    }
    return sum * dx * dy;
}

double local_abs_mass(const Field& field, double dx, double dy) {
    double sum = 0.0;
    for (int i = 1; i <= field.ny(); ++i) {
        const double* row = field.row(i);
        for (int j = 1; j <= field.nx(); ++j) sum += std::abs(row[j]);
    }
    return sum * dx * dy;
}

double cosine_decay_continuous(const Config& c, double x, double y, double t) {
    const double rate = c.kappa * pi * pi * (1 / (c.lx * c.lx) + 1 / (c.ly * c.ly));
    return c.amplitude * std::cos(pi * x / c.lx) * std::cos(pi * y / c.ly) * std::exp(-rate * t);
}

double cosine_decay_discrete_factor(const Config& c, double dt, int steps) {
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;
    const double sx = std::sin(pi * dx / (2 * c.lx));
    const double sy = std::sin(pi * dy / (2 * c.ly));
    const double lambda = 4 * c.kappa / (dx * dx) * sx * sx + 4 * c.kappa / (dy * dy) * sy * sy;
    return std::pow(1 - lambda * dt, steps);
}

CosineErrors cosine_decay_errors(const Field& field, const Config& c, const Block& b, double dt, int steps) {
    const double dx = c.lx / c.nx;
    const double dy = c.ly / c.ny;
    const double factor = cosine_decay_discrete_factor(c, dt, steps);
    const double t = steps * dt;
    const double scale = c.amplitude != 0.0 ? std::abs(c.amplitude) : 1.0;

    CosineErrors e;
    for (int i = 1; i <= b.ny; ++i) {
        const double y = cell_centre_y(b.row_offset + i, dy, c.ly);
        for (int j = 1; j <= b.nx; ++j) {
            const double x = cell_centre_x(b.col_offset + j, dx);
            const double value = field(i, j);
            e.discrete = std::max(e.discrete, std::abs(value - factor * initial_value(c, x, y)));
            e.continuous = std::max(e.continuous, std::abs(value - cosine_decay_continuous(c, x, y, t)));
        }
    }
    e.discrete /= scale;
    e.continuous /= scale;
    return e;
}
