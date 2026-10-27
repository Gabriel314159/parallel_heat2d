// Unit test for src/common/grid.hpp. Built and run by CTest: ctest --test-dir build

#include "common/grid.hpp"

#include <iostream>
#include <set>
#include <string>
#include <utility>  // std::swap

namespace {

int failures = 0;

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

}  // namespace

int main() {
    const int ny = 3;
    const int nx = 4;

    std::cout << "Field\n";
    {
        Field c(ny, nx);
        check(c.values().size() == static_cast<std::size_t>((ny + 2 * halo) * (nx + 2 * halo)),
              "size = (ny+2)(nx+2) = 30 doubles");
        check(c.stride() == nx + 2 * halo, "stride = nx+2 = 6");

        bool all_zero = true;
        for (double v : c.values()) all_zero = all_zero && v == 0.0;
        check(all_zero, "starts at 0.0");

        // Write a different number into every cell, ghosts included.
        for (int i = 1 - halo; i <= ny + halo; ++i)
            for (int j = 1 - halo; j <= nx + halo; ++j)
                c(i, j) = 10.0 * i + j;

        check(c(0, 0) == 0.0 && c(2, 3) == 23.0 && c(ny + 1, nx + 1) == 45.0,
              "read back what was written");

        // Row-major: the vector must be 0,1,2,..,5, 10,11,..,15, 20,..., 40,..,45.
        bool row_major = true;
        int k = 0;
        for (int i = 1 - halo; i <= ny + halo; ++i)
            for (int j = 1 - halo; j <= nx + halo; ++j)
                row_major = row_major && c.values()[k++] == 10.0 * i + j;
        check(row_major, "stored row by row, top ghost row first");

        // A const reference can only use the const member functions.
        const Field& view = c;
        check(view(1, 1) == 11.0 && view.ny() == ny && view.nx() == nx, "read through a const Field&");

        c.fill(7.0);
        check(c(1, 1) == 7.0 && c(0, nx + 1) == 7.0, "fill sets every value");
    }

    std::cout << "Swap\n";
    {
        Field current(ny, nx);
        Field next(ny, nx);
        current.fill(1.0);
        next.fill(2.0);
        const double* current_storage = current.values().data();
        // std::swap exchanges the two vectors' internal pointers: no values are
        // copied, so swapping a 4096 x 4096 field is as cheap as this one.
        std::swap(current, next);
        check(current(1, 1) == 2.0 && next(1, 1) == 1.0, "std::swap exchanges the fields");
        check(next.values().data() == current_storage, "swap moves storage instead of copying it");
    }

    std::cout << "FaceVelocities\n";
    {
        FaceVelocities u(ny, nx);
        std::set<double*> east_slots;
        std::set<double*> north_slots;
        for (int i = 1; i <= ny; ++i)
            for (int j = 0; j <= nx; ++j) east_slots.insert(&u.east(i, j));
        for (int i = 1; i <= ny + 1; ++i)
            for (int j = 1; j <= nx; ++j) north_slots.insert(&u.north(i, j));
        check(east_slots.size() == static_cast<std::size_t>(ny * (nx + 1)),
              "every east face has its own slot (15)");
        check(north_slots.size() == static_cast<std::size_t>((ny + 1) * nx),
              "every north face has its own slot (16)");

        u.east(2, 0) = 0.5;
        u.north(ny + 1, nx) = -0.25;
        const FaceVelocities& view = u;
        check(view.east(2, 0) == 0.5 && view.north(ny + 1, nx) == -0.25, "write and read faces");
        check(&u.east(1, 1) == &u.east(1, 0) + 1, "east faces stored row by row");
        check(&u.north(2, 1) == &u.north(1, 1) + nx, "north faces stored row by row");
    }

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
