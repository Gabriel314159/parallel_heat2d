// Storage for the 3D pigment field (docs/design_3d.md, "Model A: 3D pigment
// in the 2D flow", "Grid").
//
// The 2D Field with a third index: (k, i, j) = (layer, row, column). Layers
// count down from the surface (k = 1 is the top layer, at z = lz), as rows
// count down from the north wall. One ghost layer surrounds the block on all
// six sides. Storage is one contiguous array, layer by layer, row by row,
// with x (the column j) varying fastest:
//
//   index(k, i, j) = ((k - 1 + halo) * (ny + 2 halo) + (i - 1 + halo)) * (nx + 2 halo) + (j - 1 + halo)
//
// so row(k, i)[j] walks along x exactly as Field::row(i)[j] does in 2D.

#pragma once

#include <cassert>
#include <cstddef>

#include "common/grid.hpp"

class Field3D {
public:
    Field3D(int nz, int ny, int nx)
        : nz_(nz), ny_(ny), nx_(nx),
          values_(static_cast<std::size_t>(nz + 2 * halo) * static_cast<std::size_t>(ny + 2 * halo) *
                  static_cast<std::size_t>(nx + 2 * halo)) {
        // First touch row by row, as the update loops share out the rows of all layers.
        zero_rows(values_.data(), (nz_ + 2 * halo) * (ny_ + 2 * halo), static_cast<std::size_t>(stride()));
    }

    int nz() const { return nz_; }
    int ny() const { return ny_; }
    int nx() const { return nx_; }

    // Doubles in one stored row, and in one stored layer, ghosts included.
    int stride() const { return nx_ + 2 * halo; }
    std::size_t layer_stride() const {
        return static_cast<std::size_t>(ny_ + 2 * halo) * static_cast<std::size_t>(stride());
    }

    double& operator()(int k, int i, int j) { return values_[index(k, i, j)]; }
    double operator()(int k, int i, int j) const { return values_[index(k, i, j)]; }

    // Pointer to cell (k, i, 0), so that row(k, i)[j] is c(k, i, j).
    double* row(int k, int i) { return &values_[index(k, i, 0)]; }
    const double* row(int k, int i) const { return &values_[index(k, i, 0)]; }

    void fill(double value) { values_.assign(values_.size(), value); }

    const Storage& values() const { return values_; }
    double* data() { return values_.data(); }
    const double* data() const { return values_.data(); }
    std::size_t size() const { return values_.size(); }

private:
    std::size_t index(int k, int i, int j) const {
        assert(k >= 1 - halo && k <= nz_ + halo);
        assert(i >= 1 - halo && i <= ny_ + halo);
        assert(j >= 1 - halo && j <= nx_ + halo);
        return (static_cast<std::size_t>(k - 1 + halo) * static_cast<std::size_t>(ny_ + 2 * halo) +
                static_cast<std::size_t>(i - 1 + halo)) *
                   static_cast<std::size_t>(stride()) +
               static_cast<std::size_t>(j - 1 + halo);
    }

    int nz_;
    int ny_;
    int nx_;
    Storage values_;
};
