// Initial conditions, the no-flux wall rule, mass diagnostics and the
// analytic cosine_decay reference. See docs/design.md, sections
// "Boundary conditions", "Initial conditions" and "Diagnostics".

#pragma once

#include "common/config.hpp"
#include "common/grid.hpp"

// The whole global grid as one block, with walls on all four sides.
Block whole_domain(const Config& c);

// Cell-centre coordinates of global cell (ig, jg); rows are numbered top-down.
inline double cell_centre_x(int jg, double dx) { return (jg - 0.5) * dx; }
inline double cell_centre_y(int ig, double dy, double ly) { return ly - (ig - 0.5) * dy; }

// The initial value of the configured case at the physical point (x, y).
double initial_value(const Config& c, double x, double y);

// Set every interior cell of a local block from its global cell-centre
// coordinates, so that every decomposition starts from the same field.
void set_initial_condition(Field& field, const Config& c, const Block& b);

// No-flux walls: on every side of the block that lies on a physical wall,
// copy the interior cells into the ghost cells mirrored across the wall
// (c_ghost = c_adjacent). Call at the start of every step.
void fill_wall_ghosts(Field& field, const Block& b);

// Interior sums of c * dx * dy and |c| * dx * dy for one block. The parallel
// solver adds these up across ranks.
double local_mass(const Field& field, double dx, double dy);
double local_abs_mass(const Field& field, double dx, double dy);

// cosine_decay, docs/design.md "Analytical reference: cosine decay".
// The continuous solution at time t:
double cosine_decay_continuous(const Config& c, double x, double y, double t);
// The exact discrete solution is (1 - lambda dt)^steps times the initial field:
double cosine_decay_discrete_factor(const Config& c, double dt, int steps);

// Largest |c - reference| over one block's interior, divided by |A0|, against
// the exact discrete and the continuous solution after `steps` steps.
struct CosineErrors {
    double discrete = 0.0;
    double continuous = 0.0;
};
CosineErrors cosine_decay_errors(const Field& field, const Config& c, const Block& b, double dt, int steps);
