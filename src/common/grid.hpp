// Storage for the scalar field and the face velocities.
// See docs/design.md, "Grid and indexing" and "Data storage and shared computation".
//
// This header has no .cpp: every function is short and is called for every
// cell at every step, so it is defined right here, where the compiler can
// inline it into the loops. Functions defined inside a class body are
// automatically "inline".

#pragma once

#include <algorithm>  // std::fill_n
#include <cassert>    // assert: checks that are compiled out in Release builds
#include <cstddef>    // std::size_t
#include <memory>     // std::allocator
#include <new>        // placement new
#include <utility>    // std::forward
#include <vector>

// Width of the ghost/halo layer on every side. The upwind + five-point
// stencil only needs 1. Use this name instead of a literal 1 everywhere.
constexpr int halo = 1;

// Allocator that leaves new elements unset. std::vector would write the zeros
// from one thread, and Linux puts a page on the socket of the thread that
// writes it first; zero_rows does the first write instead.
template <class T>
struct UnsetAllocator {
    using value_type = T;
    UnsetAllocator() = default;
    template <class U>
    UnsetAllocator(const UnsetAllocator<U>&) noexcept {}
    T* allocate(std::size_t n) { return std::allocator<T>().allocate(n); }
    void deallocate(T* p, std::size_t n) noexcept { std::allocator<T>().deallocate(p, n); }
    template <class U>
    void construct(U* p) { ::new (static_cast<void*>(p)) U; }
    template <class U, class... Args>
    void construct(U* p, Args&&... args) { ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...); }
    bool operator==(const UnsetAllocator&) const noexcept { return true; }
    bool operator!=(const UnsetAllocator&) const noexcept { return false; }
};
using Storage = std::vector<double, UnsetAllocator<double>>;

// Zero `rows` rows of `width` doubles, shared out like the solver loops
// (OpenMP static schedule): each thread first touches the rows it computes on.
inline void zero_rows(double* data, int rows, std::size_t width) {
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int r = 0; r < rows; ++r) std::fill_n(data + static_cast<std::size_t>(r) * width, width, 0.0);
}

// Where a local block of cells sits in the global grid, and which of its four
// sides lie on a physical wall. The serial solver uses one block covering the
// whole grid with walls on all sides; the parallel solver has one per rank.
struct Block {
    int ny = 0;          // local rows
    int nx = 0;          // local columns
    int row_offset = 0;  // local row i is global row row_offset + i
    int col_offset = 0;  // local column j is global column col_offset + j
    bool wall_top = true;
    bool wall_bottom = true;
    bool wall_left = true;
    bool wall_right = true;
};

// A 2D field of doubles with a ghost layer, stored row by row in ONE
// contiguous std::vector (row-major order).
//
// Indices follow docs/design.md: rows i = 1..ny and columns j = 1..nx are the
// interior cells; i = 0, i = ny+1, j = 0, j = nx+1 are the ghost cells
// (for halo = 1). In general i runs from 1-halo to ny+halo.
//
// A "class" is a struct whose members are private by default: code outside the
// class can only use the public member functions. That lets Field guarantee
// that the vector always has the right size.
class Field {
public:
    // Constructor: runs when a Field is created, e.g.  Field c(ny, nx);
    // The part after ":" initializes the members before the body runs.
    Field(int ny, int nx)
        : ny_(ny), nx_(nx),
          values_(static_cast<std::size_t>(ny + 2 * halo) * static_cast<std::size_t>(nx + 2 * halo)) {
        zero_rows(values_.data(), ny_ + 2 * halo, static_cast<std::size_t>(stride()));
    }

    // "const" after the parentheses: calling this does not modify the Field.
    int ny() const { return ny_; }
    int nx() const { return nx_; }

    // Number of doubles in one stored row, ghosts included.
    int stride() const { return nx_ + 2 * halo; }

    // Element access with round brackets:  c(i, j) = 1.0;  double x = c(i, j);
    // Returning "double&" (a reference) is what makes  c(i, j) = 1.0  write
    // into the vector instead of into a temporary copy.
    double& operator()(int i, int j) { return values_[index(i, j)]; }

    // The read-only version, used when the Field itself is const.
    double operator()(int i, int j) const { return values_[index(i, j)]; }

    // Pointer to cell (i, 0), so that  row(i)[j]  is c(i, j) for every stored j.
    // The solver loops use these to walk along a row without recomputing the
    // index of every cell.
    double* row(int i) { return &values_[index(i, 0)]; }
    const double* row(int i) const { return &values_[index(i, 0)]; }

    void fill(double value) { values_.assign(values_.size(), value); }

    // Raw access to the whole vector, ghosts included. Used for output and tests.
    const Storage& values() const { return values_; }

    // The storage as one array, ghosts included: (i, j) is at i * stride() + j.
    double* data() { return values_.data(); }
    const double* data() const { return values_.data(); }
    std::size_t size() const { return values_.size(); }

private:
    // Position of cell (i, j) in values_. The first stored row is i = 1-halo
    // and the first stored column is j = 1-halo.
    std::size_t index(int i, int j) const {
        assert(i >= 1 - halo && i <= ny_ + halo);
        assert(j >= 1 - halo && j <= nx_ + halo);
        return static_cast<std::size_t>(i - 1 + halo) * static_cast<std::size_t>(stride()) +
               static_cast<std::size_t>(j - 1 + halo);
    }

    // Trailing underscore: a common convention for private member names.
    int ny_;
    int nx_;
    Storage values_;
};

// Normal velocities on the cell faces of an ny x nx block (no ghosts needed).
//   east(i, j): east face of cell (i, j),  i = 1..ny,    j = 0..nx
//               (j = 0 is the west edge of the block, j = nx its east edge)
//   north(i, j): north face of cell (i, j), i = 1..ny+1, j = 1..nx
//               (i = 1 is the top edge of the block, i = ny+1 its bottom edge,
//                which is also "the north face of the row below the block")
// So the west face of (i, j) is east(i, j-1) and its south face is north(i+1, j).
class FaceVelocities {
public:
    FaceVelocities(int ny, int nx)
        : ny_(ny),
          nx_(nx),
          east_(static_cast<std::size_t>(ny) * static_cast<std::size_t>(nx + 1)),
          north_(static_cast<std::size_t>(ny + 1) * static_cast<std::size_t>(nx)) {
        zero_rows(east_.data(), ny, static_cast<std::size_t>(nx + 1));
        zero_rows(north_.data(), ny + 1, static_cast<std::size_t>(nx));
    }

    int ny() const { return ny_; }
    int nx() const { return nx_; }

    double& east(int i, int j) { return east_[east_index(i, j)]; }
    double east(int i, int j) const { return east_[east_index(i, j)]; }
    double& north(int i, int j) { return north_[north_index(i, j)]; }
    double north(int i, int j) const { return north_[north_index(i, j)]; }

    // Row pointers for the solver loops:
    //   east_row(i)[j]       is east(i, j),  j = 0..nx
    //   north_row(i)[j - 1]  is north(i, j), j = 1..nx
    double* east_row(int i) { return &east_[east_index(i, 0)]; }
    const double* east_row(int i) const { return &east_[east_index(i, 0)]; }
    double* north_row(int i) { return &north_[north_index(i, 1)]; }
    const double* north_row(int i) const { return &north_[north_index(i, 1)]; }

private:
    std::size_t east_index(int i, int j) const {
        assert(i >= 1 && i <= ny_);
        assert(j >= 0 && j <= nx_);
        return static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(nx_ + 1) + static_cast<std::size_t>(j);
    }

    std::size_t north_index(int i, int j) const {
        assert(i >= 1 && i <= ny_ + 1);
        assert(j >= 1 && j <= nx_);
        return static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(nx_) + static_cast<std::size_t>(j - 1);
    }

    int ny_;
    int nx_;
    Storage east_;
    Storage north_;
};
