# Technical references

These sources informed the planned design. Record any additional algorithmic or implementation references used during development.

- ENCCS, [Heat diffusion mini-app](https://enccs.github.io/openmp-gpu/miniapp/): five-point finite differences, explicit time stepping and field storage.
- ENCCS, [Non-blocking point-to-point](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/): halo exchange, nonblocking messages and communication/computation overlap.
- ENCCS, [MPI and threads in practice](https://enccs.github.io/intermediate-mpi/mpi-and-threads-pt2/): hybrid MPI/OpenMP organization and thread placement.
- MPI Forum, [MPI 4.1 standard](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/mpi41-report.htm): message passing, thread support, topology and timing semantics.
- CMake, [FindMPI](https://cmake.org/cmake/help/latest/module/FindMPI.html) and [FindOpenMP](https://cmake.org/cmake/help/latest/module/FindOpenMP.html): imported build targets and dependency discovery.

The source code will use standard MPI C calls from C++, not the removed legacy MPI C++ bindings. Platform module names, launch commands and queue settings must come from current Atlas/Vanda documentation.
