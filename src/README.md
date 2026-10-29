# C++ sources

The shared components (`common/`), the serial solver (`serial/`) and the MPI + OpenMP solver (`parallel/`) are implemented and tested. The model, discretization and boundary rule are specified in `docs/design.md`. Unit tests are in `tests/` and run with `ctest --test-dir build`.

| Directory | Contents |
| --- | --- |
| `common/` | Components shared by both solvers: configuration, grid, the flow (vorticity equation and Poisson solver), physics, stencil and output. Independent of MPI and OpenMP. |
| `serial/` | The serial driver and solver, the reference for correctness. |
| `parallel/` | The MPI + OpenMP driver, solver, decomposition and halo exchange. The solver supplies the global sums and halo exchanges of the shared CG Poisson solver. |

Source includes can use paths such as `common/grid.hpp`; CMake exposes the `src` directory to all targets.

## Shared components (`common/`)

| File | Responsibility |
| --- | --- |
| `config.hpp`, `config.cpp` | CLI parsing, parameter validation and the fixed time step from `cfl` and `--max-speed`. Also the output options (`--output`, snapshots) and the `--help` text. |
| `stability.hpp`, `stability.cpp` | The three stability limits (pigment, vorticity advection, viscosity) of `docs/design.md`, "Time step": the largest stable `dt` for given speeds, the check of measured speeds during the run, and the largest face speeds of a block. |
| `grid.hpp` | `Field`: a contiguous `std::vector<double>`, local dimensions, a halo of width `halo` on every side, `(i, j)` access and row pointers for the loops. `FaceVelocities`: the two face-velocity arrays. `Block`: where a local block sits in the global grid and which of its sides are walls. |
| `corners.hpp` | Corner fields ($\psi$, $\zeta$) stored in a `Field`: corner positions, the range of unknown corners of a block, and zeroing the ghost corners on the top and left walls. |
| `flow.hpp`, `flow.cpp` | One row of the vorticity right-hand side $R = J(\psi,\zeta) + F + \nu\nabla^2\zeta$ (Arakawa Jacobian), the Adams–Bashforth weights and update, the wind profile and strength, the face velocities from $\psi$, and the energy and enstrophy sums. |
| `poisson.hpp` | The conjugate gradient Poisson solver, matrix-free and warm-started. Header-only, so the parallel solver compiles its loops with OpenMP; the global sum and the halo exchange come in through a small `PoissonComm` interface (`SerialComm` for one process). |
| `velocity.hpp`, `velocity.cpp` | An analytic streamfunction (the former prescribed double gyre) and its per-face velocities, kept as a reference for the face-velocity tests. |
| `physics.hpp`, `physics.cpp` | Initial conditions (`pigment_drop` Gaussian, `constant`, `cosine_decay`), the no-flux wall ghost fill, the analytic `cosine_decay` reference (continuous and exact-discrete), and the local mass sum $\Delta x\,\Delta y\sum c_{i,j}$. |
| `stencil.hpp` | An inline upwind flux update of one pigment cell from its five cell values and its four face velocities; no communication or threading inside this function. |
| `io.hpp`, `io.cpp` | Final fields as CSV or `.npy`, snapshot folders (`meta.json`, `index.csv`, pigment, $\psi$ and vorticity as `.npy`, optionally coarsened), and the run summary: a short report plus one `RESULT {json}` line with the parameters, `dt`, timings, mass error, energy, enstrophy and measured speeds. |

### Stencil

For cell $(i,j)$, with face velocities $U_e(i,j)$, $U_e(i,j-1)$, $V_n(i,j)$, $V_n(i+1,j)$:

$$
c^{\,n+1}_{i,j} = c^{\,n}_{i,j}
- \frac{\Delta t}{\Delta x}\big(F_e(i,j) - F_e(i,j-1)\big)
- \frac{\Delta t}{\Delta y}\big(G_n(i,j) - G_n(i+1,j)\big),
$$

where each face flux is the upwind advective flux minus the central diffusive flux, exactly as in `docs/design.md`. One helper, `face_flux(c_minus, c_plus, vel, kappa_h)`, computes the flux through a face from the cell on its minus side (west or south) to the cell on its plus side (east or north), and both cells that share the face call it with the same arguments in the same order. Then interior fluxes cancel bit-for-bit and the mass is conserved to rounding, in serial and across ranks.

With no wind the flow stays at rest, all face velocities are zero, and the update reduces to five-point diffusion.

### The flow

Each solver owns the corner fields $\psi$ and $\zeta$ and three right-hand-side fields $R^n$, $R^{n-1}$, $R^{n-2}$, all of the pigment field's shape (see `corners.hpp`). One flow step (`docs/design.md`, "One time step"):

1. Exchange $\zeta$'s halos, including the corner ghosts (parallel only; in serial the wall ghosts simply stay 0).
2. Compute $R^n$ at the unknown corners and advance $\zeta$ with the Adams–Bashforth weights of this step; then rotate the three $R$ fields by swapping, without copying.
3. Solve the Poisson equation for $\psi$ with CG, starting from the previous $\psi$ (`PoissonCG`; in parallel with `MPI_Allreduce` sums and the halo exchange).
4. Exchange $\psi$'s halos with corners (parallel only) and compute all face velocities from it.

Face velocities are stored in arrays, as for the pigment stencil, and refreshed once per step. Faces on a wall come out exactly 0, because $\psi = 0$ at both of their corners.

The time step is fixed for the run. Every `check_interval` steps the largest face speeds are measured (one `MPI_Allreduce` in parallel) and checked against the stability limits; a run whose speeds outgrow `--max-speed` stops with a message.

### Invariants

- Grid access uses `(row, column)`; `nx` counts columns and `ny` counts rows. Rows are numbered top-down.
- The grid is cell-centred: $\Delta x = l_x/n_x$, $\Delta y = l_y/n_y$, and cell $(i,j)$ has centre $\big((j-\tfrac12)\Delta x,\ l_y-(i-\tfrac12)\Delta y\big)$ in global indices.
- The halo width is a named constant `halo` (currently 1), not a literal `1` spread through the code. Allocate `(local_ny+2*halo)*(local_nx+2*halo)` doubles per field. A higher-order advection scheme would then only change `halo` and the stencil.
- Keep two equally sized fields and swap their underlying buffers after each step.
- Physical-wall ghosts depend on the current interior, so they are **refilled at the start of every step** by the boundary rule, not set once at initialization. Process halos are filled by communication.
- Corner ghost cells are never read by the five-point update and need no value.
- Initialize each local cell from its global cell-centre coordinates so all decompositions represent the same problem.
- Keep shared components independent of MPI and OpenMP; parallel-specific gathering and reductions belong in the parallel layer.
- Do not put file output inside the point-update kernel.
- Keep `double` arithmetic and the same expression order across baseline and parallel updates.

## Serial solver (`serial/`)

Files: `main.cpp`, `solver.hpp`, `solver.cpp`. Run `heat_serial --help` for the options.

The driver parses parameters, allocates the full domain with its ghost layer, the face-velocity arrays and the flow fields, spins the flow up from rest, releases the pigment, runs a fixed number of time steps, optionally writes the final field and snapshots, and reports the timings, the mass error and the flow diagnostics.

During spin-up each step is a flow step only. After the release, each time step $n$ is:

1. Refill the physical-wall ghosts of the pigment from its interior (no-flux: ghost = adjacent interior cell).
2. Update every pigment cell with the shared stencil and the face velocities of $\psi^n$, then swap the current and next fields.
3. Advance the flow to $\zeta^{n+1}$, $\psi^{n+1}$ and the new face velocities (see "The flow" above).
4. Every `check_interval` steps, check the measured speeds against the stability limits.

The solver must not depend on MPI or OpenMP. Use the same problem, time step, number of steps and precision as the parallel versions.

Timing covers the whole time-stepping loop, including the ghost fill and the flow step; the spin-up and the snapshot writing are measured separately, and the time spent in the Poisson solve is reported as well. Allocation, initialization, file output and validation are separate phases. The total mass $M$ is computed before and after the loop, outside the timed region; with no-flux walls its relative change (`mass_error`, see `docs/experiment_plan.md`) stays at rounding level, about $10^{-13}$ after tens of thousands of steps. Do not infer correctness from a heatmap alone.

## Parallel solver (`parallel/`)

Run `heat_parallel --help` for the options. On every correctness case it agrees with the serial field to about $10^{-13}$ (the CG sums add in a different order); `ctest` checks this (`correctness_parallel`).

| File | Responsibility |
| --- | --- |
| `main.cpp` | MPI initialization (`MPI_THREAD_FUNNELED`, `MPI_ERRORS_RETURN`), parameter handling, solver invocation and finalization. Any exception prints the rank and ends in `MPI_Abort`. |
| `solver.hpp`, `solver.cpp` | The thread count (`configure_threads`), spin-up and time loop, physical-wall ghost fill, the flow step, OpenMP worksharing, strategy selection, the speed check, timing and the global sums. |
| `decomposition.hpp`, `decomposition.cpp` | Local sizes and global offsets, remainder distribution, process topology (Cartesian communicator without reordering) and physical-edge ownership; the node count; `mpi_check`, which turns an MPI error code into an exception. |
| `halo.hpp`, `halo.cpp` | Blocking exchange and start/finish operations for nonblocking exchange of the pigment, and the blocking exchange with corners for $\psi$ and $\zeta$. |

### Implementation

1. Request `MPI_THREAD_FUNNELED` and verify the returned support level. Only the initializing thread calls MPI. Set `MPI_ERRORS_RETURN` on the world and the Cartesian communicator and wrap every MPI call in `mpi_check`; `main` turns any exception into a message with the rank and `MPI_Abort`. Choose the threads per rank (`--threads`, else `OMP_NUM_THREADS`, else 1), check that all ranks agree, and count the nodes (`MPI_Comm_split_type` with `MPI_COMM_TYPE_SHARED`).
2. Split global rows across ranks, including any remainder. Then extend the same interface to a 2D Cartesian process grid. The decomposition is unchanged by the advection–diffusion model.
3. Store local fields with a halo of width `halo` on all sides, plus the local face-velocity arrays.
4. Every step, fill each side of the local block in one of two ways:
   - **Internal edge** (a neighbouring rank exists): MPI halo exchange.
   - **Physical edge** (the side lies on a global wall): local ghost fill by the boundary rule, e.g. the no-flux mirror fill, ghost = adjacent interior cell. This is refilled every step, because it depends on the current interior.
5. Implement blocking halo exchange using a deadlock-safe communication pattern such as `MPI_Sendrecv`.
6. After the pigment update, advance the flow: exchange $\zeta$ with corners, compute the right-hand side and the Adams–Bashforth update on the local unknown corners, solve the Poisson equation with CG (one halo exchange and two `MPI_Allreduce` sums per iteration), exchange $\psi$ with corners, and compute the local face velocities. Both ranks next to a shared face compute it from the same $\psi$ values, so they agree exactly.
7. Parallelize the pigment update, the vorticity rows, the face velocities and the CG loops with OpenMP, with one writer per destination value.
8. Overlap (pigment only): post receives and sends, fill the physical-wall ghosts, update the strict local interior, wait for all requests, then update the remaining local boundary cells exactly once. The flow exchanges are blocking.
9. Complete outstanding sends and receives before swapping or reusing their buffers. Do not update a receive region while communication is pending.

### Strict interior under overlap

The strict interior is the set of owned cells whose stencil reads no halo cell that is still being received. With halo width $h$ (currently 1), that means cells at least $h$ cells away from every **internal** edge. Cells next to a **physical** edge may be in the strict interior, because their ghosts are filled locally before the update. If `halo` is increased later, this region shrinks accordingly; derive it from `halo`, do not hard-code an offset of 1.

At one or two cells per local dimension, the overlap interior may be empty. Use a correct boundary-region partition without duplicate updates or out-of-range accesses. Reject decompositions that assign zero-sized local domains.

### Communication details

Rows are contiguous: they are sent straight from the field with an MPI vector datatype (`halo` rows of `nx` interior values). Columns are packed into contiguous buffers before sending and unpacked after receiving. Blocking and overlapping exchange use the same scheme, so comparing them changes only the timing of the messages. For the pigment, only interior values are sent and only edge ghost cells are written. For $\psi$ and $\zeta$, the corner exchange first sends the columns and then whole rows including their ghost columns, which carries the diagonal neighbours' values into the corner ghosts.

The Cartesian communicator is created with `reorder = 0`, so its ranks are the world ranks and the launcher's placement of ranks on cores is kept. Still use the Cartesian communicator for everything after the layout is made (neighbors, gathers, reductions), and determine neighbors and gather metadata with its own rank numbering.

### Diagnostics and timing

The total mass $M = \Delta x\,\Delta y\sum c_{i,j}$, the energy and the enstrophy need a local sum followed by an `MPI_Allreduce`. They are computed before and after the timed loop and at snapshots, never every step inside the timed region, and differ from the serial values in the last digits because the summation order differs. The speed check every `check_interval` steps needs one `MPI_Allreduce` with `MPI_MAX`, which is exact.

Gather only small fields for validation or selected visual output, outside solver timing. Use a barrier before the timed region and reduce local elapsed times with `MPI_MAX` after it. The timed loop includes the halo exchanges, the physical-wall ghost fill and the whole flow step. A global barrier is not required after every time step when neighbor exchanges already enforce the required dependencies.
