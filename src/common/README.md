# Shared C++ components to implement

This directory intentionally contains no solver source files yet.

| Planned file | Responsibility |
| --- | --- |
| `config.hpp`, `config.cpp` | CLI parsing, parameter validation and global time-step calculation. |
| `grid.hpp` | A contiguous `std::vector<double>`, local dimensions, one-cell halos and a row/column indexing API. |
| `physics.hpp`, `physics.cpp` | Initial conditions, physical boundary assignment and analytical reference fields. |
| `stencil.hpp` | An inline update of one point using `rx` and `ry`; no communication or threading inside this function. |
| `io.hpp`, `io.cpp` | Serial matrix export, metadata and benchmark record formatting. |

## Invariants

- Grid access uses `(row, column)`; `nx` counts columns and `ny` counts rows.
- Allocate `(local_ny+2)*(local_nx+2)` doubles per field, including halos.
- Keep two equally sized fields and swap their underlying buffers after each step.
- Initialize physical boundaries in both buffers. Process halos are filled by communication.
- Initialize each local point from its global coordinates so all decompositions represent the same problem.
- Keep shared components independent of MPI and OpenMP; parallel-specific gathering belongs in the parallel layer.
- Do not put file output inside the point-update kernel.
- Keep `double` arithmetic and the same expression order across baseline and parallel updates where practical.

Source includes can use paths such as `common/grid.hpp`; CMake exposes the `src` directory to all targets.
