# Serial solver

Files: `main.cpp`, `solver.hpp`, `solver.cpp`.

The driver parses parameters, allocate the full plate, initialize both buffers, run a fixed number of time steps, optionally write the final field, and report solver time.

The solver visits all interior rows and columns in a deterministic loop and call the shared stencil update. It must not depend on MPI or OpenMP. Use the same problem, time step, number of steps and precision as the parallel versions.

Timing will cover only the time-stepping loop. Allocation, initialization, file output and validation are separate phases. Do not infer correctness from a heatmap alone.
