# C++ sources

The shared components (`common/`), the serial solver (`serial/`) and the MPI + OpenMP solver (`parallel/`) are implemented and tested. The model, discretization and boundary rule are specified in `docs/design.md`. Unit tests are in `tests/` and run with `ctest --test-dir build`.

| Directory | Contents |
| --- | --- |
| `common/` | Components shared by both solvers: configuration, grid, velocity, physics, stencil and output. Independent of MPI and OpenMP. |
| `serial/` | The serial driver and solver, the reference for correctness. |
| `parallel/` | The MPI + OpenMP driver, solver, decomposition and halo exchange. |

Source includes can use paths such as `common/grid.hpp`; CMake exposes the `src` directory to all targets.

## Shared components (`common/`)

| File | Responsibility |
| --- | --- |
| `config.hpp`, `config.cpp` | CLI parsing, parameter validation and the global time step. Validates $\kappa > 0$, $\varepsilon \ge 0$, flow period $T > 0$, blob width $\sigma > 0$ and the case requirements. Computes $u_{\max}$, $v_{\max}$ from `velocity.hpp`, then $a$, $d$, and either $\Delta t = \mathrm{cfl}/(a+d)$ or checks an explicit `dt` against $\Delta t\,(a+d) \le 1$. Also computes $\mathrm{Pe}$ and $\mathrm{Pe}_{\text{grid}}$ for the run metadata, and holds the output options (`--output`, snapshots) and the `--help` text. |
| `grid.hpp` | `Field`: a contiguous `std::vector<double>`, local dimensions, a halo of width `halo` on every side, `(i, j)` access and row pointers for the loops. `FaceVelocities`: the two face-velocity arrays (see below). `Block`: where a local block sits in the global grid and which of its sides are walls. |
| `velocity.hpp`, `velocity.cpp` | The double-gyre streamfunction $\psi(x,y,t)$, the analytic speed bounds $u_{\max}$, $v_{\max}$, the per-face functions `face_velocity_east/north`, and `FlowField`, which fills a whole block's face velocities at time $t$ from global coordinates. |
| `physics.hpp`, `physics.cpp` | Initial conditions (`pigment_drop` Gaussian, `constant`, `cosine_decay`), the no-flux wall ghost fill, the analytic `cosine_decay` reference (continuous and exact-discrete), and the local mass sum $\Delta x\,\Delta y\sum c_{i,j}$. |
| `stencil.hpp` | An inline upwind flux update of one cell from its five cell values and its four face velocities; no communication or threading inside this function. |
| `io.hpp`, `io.cpp` | Final fields as CSV or `.npy`, snapshot folders (`meta.json`, `index.csv`, `.npy` fields, optionally coarsened), and the run summary: a short report plus one `RESULT {json}` line with the parameters, `dt`, timing, mass error and Péclet numbers. |

### Stencil

For cell $(i,j)$, with face velocities $U_e(i,j)$, $U_e(i,j-1)$, $V_n(i,j)$, $V_n(i+1,j)$:

$$
c^{\,n+1}_{i,j} = c^{\,n}_{i,j}
- \frac{\Delta t}{\Delta x}\big(F_e(i,j) - F_e(i,j-1)\big)
- \frac{\Delta t}{\Delta y}\big(G_n(i,j) - G_n(i+1,j)\big),
$$

where each face flux is the upwind advective flux minus the central diffusive flux, exactly as in `docs/design.md`. One helper, `face_flux(c_minus, c_plus, vel, kappa_h)`, computes the flux through a face from the cell on its minus side (west or south) to the cell on its plus side (east or north), and both cells that share the face call it with the same arguments in the same order. Then interior fluxes cancel bit-for-bit and the mass is conserved to rounding, in serial and across ranks.

The stencil replaces the old heat-only update with `rx` and `ry`. With $A = 0$ all face velocities are zero and the update reduces to the old five-point diffusion.

### Velocity storage

Face velocities are stored in arrays, not recomputed inside the stencil:

- `U` has `local_ny` rows of `local_nx+1` vertical faces, `V` has `local_ny+1` rows of `local_nx` horizontal faces.
- Each interior face is shared by two cells, so storing it evaluates $\psi$ (two `sin` calls per corner) once per face instead of twice, and keeps the stencil free of trigonometry.
- If the flow is steady ($\varepsilon = 0$ or $A = 0$), fill the arrays once before the time loop. Otherwise refill them at the start of every step at the old time $t_n = n\,\Delta t$; this cost belongs to the timed loop.
- Evaluate $\psi$ from **global** corner coordinates. Two ranks that share a face then compute the same value without communication.
- `FlowField` uses that $\psi$ is separable, $\psi = A\,s_x(x,t)\,s_y(y)$ with $s_x = \sin(\pi f(X,t))$ and $s_y = \sin(\pi Y)$: it evaluates $s_x$ once per corner column per step and $s_y$ once per corner row, so a time-dependent refill costs $O(n_x + n_y)$ sines. The arithmetic is the same as in `face_velocity_east/north`, and the unit test checks that the two agree exactly.
- Faces on a physical wall are set to exactly 0. Evaluating $\psi$ there gives about $10^{-16}$ instead, because $\sin(2\pi)$ is not exactly 0 in floating point.

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

The driver parses parameters, allocates the full domain with its ghost layer and the two face-velocity arrays, initializes the interior, runs a fixed number of time steps, optionally writes the final field and snapshots, and reports solver time and the mass error.

Each time step $n$, at time $t_n = n\,\Delta t$:

1. Refill the physical-wall ghosts of the current field from its interior by the boundary rule (no-flux: ghost = adjacent interior cell).
2. If the flow is time-dependent ($\varepsilon > 0$ and $A \ne 0$), refill the face velocities at $t_n$. For a steady flow they are filled once before the loop.
3. Visit all interior rows and columns in a deterministic loop and call the shared stencil update.
4. Swap the current and next fields.

The solver must not depend on MPI or OpenMP. Use the same problem, time step, number of steps and precision as the parallel versions.

Timing covers the whole time-stepping loop, including the ghost fill and any per-step velocity evaluation; snapshot writing is measured separately and excluded. Allocation, initialization, file output and validation are separate phases. The total mass $M$ is computed before and after the loop, outside the timed region; with no-flux walls its relative change (`mass_error`, see `docs/experiment_plan.md`) stays at rounding level, about $10^{-13}$ after tens of thousands of steps. Do not infer correctness from a heatmap alone.

## Parallel solver (`parallel/`)

Run `heat_parallel --help` for the options. On every correctness case it reproduces the serial field bit for bit; `ctest` checks this (`correctness_parallel`).

| File | Responsibility |
| --- | --- |
| `main.cpp` | MPI initialization, rank-aware errors, parameter handling, solver invocation and finalization. |
| `solver.hpp`, `solver.cpp` | Time loop, physical-wall ghost fill, per-step velocity refill, OpenMP worksharing, strategy selection, timing and the mass reduction. |
| `decomposition.hpp`, `decomposition.cpp` | Local sizes and global offsets, remainder distribution, process topology and physical-edge ownership. |
| `halo.hpp`, `halo.cpp` | Blocking exchange and start/finish operations for nonblocking exchange of the scalar field. |

### Initial implementation

1. Request `MPI_THREAD_FUNNELED` and verify the returned support level. Only the initializing thread calls MPI.
2. Split global rows across ranks, including any remainder. Then extend the same interface to a 2D Cartesian process grid. The decomposition is unchanged by the advection–diffusion model.
3. Store local fields with a halo of width `halo` on all sides, plus the local face-velocity arrays.
4. Every step, fill each side of the local block in one of two ways:
   - **Internal edge** (a neighbouring rank exists): MPI halo exchange.
   - **Physical edge** (the side lies on a global wall): local ghost fill by the boundary rule, e.g. the no-flux mirror fill, ghost = adjacent interior cell. This is refilled every step, because it depends on the current interior.
5. Implement blocking halo exchange using a deadlock-safe communication pattern such as `MPI_Sendrecv`.
6. If the flow is time-dependent ($\varepsilon > 0$, the double gyre), refill the local face velocities at $t_n$ every step. Each rank evaluates $\psi$ at global corner coordinates, so faces on a rank boundary agree on both sides without communicating velocities. Only the scalar field is exchanged.
7. Parallelize local point updates, the ghost fill and the velocity refill with OpenMP, with one writer per destination value.
8. Add overlap: post receives and sends, fill physical-wall ghosts and refill velocities (both local work), update the strict local interior, wait for all requests, then update the remaining local boundary cells exactly once.
9. Complete outstanding sends and receives before swapping or reusing their buffers. Do not update a receive region while communication is pending.

### Strict interior under overlap

The strict interior is the set of owned cells whose stencil reads no halo cell that is still being received. With halo width $h$ (currently 1), that means cells at least $h$ cells away from every **internal** edge. Cells next to a **physical** edge may be in the strict interior, because their ghosts are filled locally before the update. If `halo` is increased later, this region shrinks accordingly; derive it from `halo`, do not hard-code an offset of 1.

At one or two cells per local dimension, the overlap interior may be empty. Use a correct boundary-region partition without duplicate updates or out-of-range accesses. Reject decompositions that assign zero-sized local domains.

### Communication details

Rows are contiguous: they are sent straight from the field with an MPI vector datatype (`halo` rows of `nx` interior values). Columns are packed into contiguous buffers before sending and unpacked after receiving. Blocking and overlapping exchange use the same scheme, so comparing them changes only the timing of the messages. Only interior values are sent and only edge ghost cells are written; corner ghosts are never touched.

For each MPI communicator used, determine neighbors and gather metadata with its own rank numbering, especially if Cartesian rank reordering is enabled.

### Diagnostics and timing

The total mass $M = \Delta x\,\Delta y\sum c_{i,j}$ needs a local sum followed by an `MPI_Allreduce`. Compute it before and after the timed loop, or only at snapshot intervals, never every step inside the timed region. Compare it with the serial value using a tolerance, because the summation order differs.

Gather only small fields for validation or selected visual output, outside solver timing. Use a barrier before the timed region and reduce local elapsed times with `MPI_MAX` after it. The timed loop includes the halo exchange, the physical-wall ghost fill and any per-step velocity refill. A global barrier is not required after every time step when neighbor exchanges already enforce the required dependencies.
