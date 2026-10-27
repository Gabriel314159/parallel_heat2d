# MPI + OpenMP solver

Both decomposition strategies and both communication modes are implemented.

| File | Responsibility |
| --- | --- |
| `main.cpp` | MPI initialization, rank-aware errors, parameter handling, solver invocation and finalization. |
| `solver.hpp`, `solver.cpp` | Time loop, OpenMP worksharing, strategy selection and timing. |
| `decomposition.hpp`, `decomposition.cpp` | Local sizes and global offsets, remainder distribution, process topology and physical-edge ownership. |
| `halo.hpp`, `halo.cpp` | Blocking exchange and start/finish operations for nonblocking exchange. |

## Implementation organization

1. Request `MPI_THREAD_FUNNELED` and verify the returned support level. Only the initializing thread calls MPI.
2. Split global rows across ranks, including any remainder. Then extend the same interface to a 2D Cartesian process grid.
3. Store local fields with a one-cell halo on all sides. Retain physical boundary values where no neighbor exists.
4. Implement blocking halo exchange using a deadlock-safe communication pattern such as `MPI_Sendrecv`.
5. Parallelize local point updates with OpenMP, with one writer per destination point.
6. Add overlap: post receives and sends, update the strict local interior, wait for all requests, then update the remaining local boundary points exactly once.
7. Complete outstanding sends and receives before swapping or reusing their buffers. Do not update a receive region while communication is pending.

Rows are contiguous. Columns need an explicit packing strategy or an MPI datatype. Pick and document one before comparing communication modes so the comparisons do not accidentally change unrelated factors.

For each MPI communicator used, determine neighbors and gather metadata with its own rank numbering, especially if Cartesian rank reordering is enabled.

At one or two cells per local dimension, the overlap interior may be empty. Use a correct boundary-region partition without duplicate updates or out-of-range accesses. Reject decompositions that assign zero-sized local domains.

Gather only small fields for validation or selected visual output, outside solver timing. Use a barrier before the timed region and reduce local elapsed times with `MPI_MAX` after it. A global barrier is not required after every time step when neighbor exchanges already enforce the required dependencies.
