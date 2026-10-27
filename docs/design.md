# Design notes

Status: planned design. The numerical and parallel implementations have not been written yet.

## Model

Solve two-dimensional advection–diffusion of a passive scalar `c(x,y,t)`, for example pigment concentration, on the rectangle `[0,lx] x [0,ly]`:

```text
dc/dt + div(u*c) = kappa * (d2c/dx2 + d2c/dy2)
```

- `kappa > 0` is a constant diffusivity.
- `u = (u, v)` is a **prescribed** velocity field. It is given as a function of position and time; it is not solved for. The scalar does not affect the flow.
- The velocity is derived from a streamfunction `psi(x,y,t)` with `u = dpsi/dy` and `v = -dpsi/dx`, so it is divergence-free. Positive `v` points towards increasing `y`, which is upwards.
- There is no source term.

With `u = 0` and Dirichlet walls the model reduces exactly to the heat equation with `kappa` as the thermal diffusivity. The heat-only cases (`heated_top`, `constant`, `sine_decay`) remain supported as this special case.

All implementations solve the same discrete problem.

## Grid and indexing

Use a cell-centred finite-volume grid. Each grid value is the average of `c` over one cell, and the physical walls lie on cell faces.

- `nx`: global cells along x (columns); `ny`: global cells along y (rows).
- `dx = lx/nx`, `dy = ly/ny`.
- Grid access uses `(row, column)`; top-down row ordering is used in field files.
- For global row `i=1..ny`, the cell centre is `y = ly-(i-0.5)*dy`; for column `j=1..nx`, it is `x = (j-0.5)*dx`.
- Row `i=0`, row `i=ny+1`, column `j=0` and column `j=nx+1` are ghost cells outside the walls. The top wall is the face between rows 0 and 1, at `y = ly`.
- The east face of cell `(i,j)` is at `x = j*dx`. The north face of cell `(i,j)` is at `y = ly-(i-1)*dy`.

The five-point update never reads a corner ghost cell, so corner values need no convention. This replaces the earlier node-based grid with `dx = lx/(nx+1)` and its corner rule.

## Velocity on faces

Store or evaluate one normal velocity per cell face, computed from `psi` at the face's two end points (cell corners):

```text
U_e(i,j) = ( psi(j*dx, ly-(i-1)*dy) - psi(j*dx, ly-i*dy) ) / dy          east face, positive = +x
V_n(i,j) = -( psi(j*dx, ly-(i-1)*dy) - psi((j-1)*dx, ly-(i-1)*dy) ) / dx  north face, positive = +y
```

The west and south faces of a cell are the east face of `(i,j-1)` and the north face of `(i+1,j)`. With this construction the discrete divergence of every cell is exactly zero, apart from rounding: the four corner values cancel. Evaluate `psi` from global coordinates so that every decomposition, and both ranks sharing a face, obtain the same face velocity without communication.

For a closed container the streamfunction is zero on all four walls, so every wall-face velocity is zero and no scalar crosses a wall by advection.

### Planned flows

| `flow` | Streamfunction | Speed bounds | Walls |
| --- | --- | --- | --- |
| `none` | `psi = 0` | 0 | any |
| `uniform` | `psi = U0*y - V0*x` | `|u| <= |U0|`, `|v| <= |V0|` | Dirichlet only; for short runs where the scalar stays away from the walls |
| `rotation` | `psi = -(omega/2)*((x-xc)^2 + (y-yc)^2)` | `|u| <= |omega|*ly/2`, `|v| <= |omega|*lx/2` for a centred axis | Dirichlet only; scalar must stay inside the inscribed circle |
| `cellular` | `psi = A*sin(m*pi*x/lx)*sin(n*pi*y/ly)` | `|u| <= |A|*n*pi/ly`, `|v| <= |A|*m*pi/lx` | closed container |
| `double_gyre` | `psi = A*sin(pi*f(X,t))*sin(pi*Y)` with `X = 2x/lx`, `Y = y/ly`, `f = eps*sin(w*t)*X^2 + (1-2*eps*sin(w*t))*X` | `|u| <= |A|*pi/ly`, `|v| <= 2*|A|*pi*(1+2*eps)/lx` | closed container |

`rotation` with positive `omega` turns counterclockwise. `cellular` with `m = n = 1` is a single vortex filling the box; larger `m`, `n` give a grid of counter-rotating cells. `double_gyre` is time-periodic and produces chaotic mixing; `eps = 0` makes it steady.

Each face velocity is a difference quotient of `psi`, so it is bounded by the maximum of the corresponding continuous velocity component. The bounds above are therefore valid bounds for the discrete face velocities and can be used in the time-step rule.

## Discretization

Write the update in conservative form. For each face, compute one total flux, advective plus diffusive, using first-order upwinding for advection and a central difference for diffusion:

```text
F_e(i,j) = U_e * (U_e > 0 ? c(i,j) : c(i,j+1))  -  kappa * (c(i,j+1) - c(i,j)) / dx
G_n(i,j) = V_n * (V_n > 0 ? c(i,j) : c(i-1,j))  -  kappa * (c(i-1,j) - c(i,j)) / dy

c_new(i,j) = c(i,j) - dt/dx * (F_e(i,j) - F_e(i,j-1)) - dt/dy * (G_n(i,j) - G_n(i+1,j))
```

The north neighbour of `(i,j)` is row `i-1`, because rows are numbered top-down.

Use forward Euler in time. The flow is evaluated at the old time level `t_n` for time-dependent flows.

Every face flux must be computed by the same expression, with the same operand order, in both cells that share the face, including cells on different ranks. Interior face fluxes then cancel when summed, and the total amount of scalar is conserved up to rounding.

Upwinding is stable and non-oscillatory but adds numerical diffusion of roughly `|u|*dx/2` in x and `|v|*dy/2` in y. The leading spatial error is therefore first order in `dx` and `dy`; it is small compared with the physical diffusion only when the grid Péclet number is small. Higher-order advection schemes would need a two-cell halo and are outside the initial scope.

## Boundary conditions

Each run uses one wall type for all four walls, set by `boundary_type`.

- `dirichlet`: a prescribed wall value `g` per wall (`top`, `bottom`, `left`, `right`). Before each step, set every ghost cell to `2*g - c_adjacent`, so that the average of the ghost and its adjacent interior cell equals `g` on the wall face.
- `noflux`: no scalar crosses a wall. Before each step, set every ghost cell to the value of its adjacent interior cell, `c_ghost = c_adjacent`. With zero wall-face velocity, both the advective and diffusive wall fluxes are then zero.

Ghost cells on physical walls now depend on the current interior values. They must be **refilled at the start of every step**, not only at initialization. Every implementation uses the same rule.

`noflux` requires a closed-container flow (`none`, `cellular` or `double_gyre`). `uniform` and `rotation` require `dirichlet`, with wall values equal to the background value so that inflowing fluid carries the background value.

## Initial conditions

Initialize interiors from cell-centre global coordinates so that a change in decomposition cannot change the initial condition.

| `case` | Initial field | Notes |
| --- | --- | --- |
| `heated_top` | `c0` everywhere | Heat-only demonstration: `flow = none`, Dirichlet walls with a hot top. |
| `constant` | `c0` everywhere | All Dirichlet wall values equal `c0`, or `noflux`. Valid with any flow. |
| `sine_decay` | `A*sin(pi*x/lx)*sin(pi*y/ly)` | `flow = none`, Dirichlet zero walls. Analytical reference. |
| `pigment_drop` | `c0 + A*exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))` | Default advection–diffusion case. |

## Time step

Define the advective and diffusive rates:

```text
a = umax/dx + vmax/dy                    (umax, vmax from the flow's speed bounds)
d = 2*kappa*(1/dx^2 + 1/dy^2)
```

Validate positive dimensions, positive `kappa`, positive `dt` and

```text
dt * (a + d) <= 1
```

For a discretely divergence-free flow this keeps every update weight non-negative in cells away from Dirichlet walls. Each new value is then a convex combination of old values, so no new maximum or minimum is created there and the scheme is stable. With `u = 0` the condition reduces to the earlier `rx + ry <= 0.5`, where `rx = kappa*dt/dx^2` and `ry = kappa*dt/dy^2`.

For automatic time stepping, `cfl` in `(0,1]` is a safety fraction:

```text
dt = cfl / (a + d)
```

A Dirichlet ghost value `2*g - c_adjacent` can lie outside the range of the initial and wall data. In cells next to a Dirichlet wall with zero wall-face velocity, the update weights stay non-negative only under the stricter condition `dt * (a + 3*kappa*(1/dx^2 + 1/dy^2)) <= 1`, so the maximum-principle check is asserted only when that holds or when walls are `noflux`.

Use fixed steps for primary scaling tests. Comparisons at different grid resolutions over the same physical duration must recompute stable time steps and the corresponding step counts; record the actual final time.

Record the Péclet numbers for each run: `Pe = U*L/kappa`, using a characteristic speed `U` and length `L`, and the grid Péclet number `Pe_grid = max(umax*dx, vmax*dy)/kappa`. When `Pe_grid` is well above 2, upwind numerical diffusion is comparable to or larger than the physical diffusion, and fine filaments are not resolved.

## Data storage and shared computation

Use contiguous `std::vector<double>` storage. Each local field includes a one-cell halo on every side, and each solver owns a current and next field. The row stride is `local_nx+2`.

For steady flows, store local face velocities: `U` has `local_ny` rows of `local_nx+1` east/west faces and `V` has `local_ny+1` rows of `local_nx` north/south faces, all computed once from global coordinates. For time-dependent flows, recompute the local face velocities at the start of each step; this cost is part of the measured solver loop. Record which approach is used, because it changes the memory traffic and arithmetic per cell.

The shared stencil computes a single destination value from the five cell values and the four face velocities of that cell. It contains no MPI calls, OpenMP regions, allocation or output. The serial and parallel solvers control the loops. This separation makes formula changes explicit and keeps the serial executable free from MPI/OpenMP dependencies.

Initialize interiors by global coordinates. Fill physical-wall ghosts by the boundary rule and process halos by exchange before their values are used in each step.

## Domain decomposition

### Rows

Distribute `ny` rows over the ranks. With `base=ny/p` and `extra=ny%p`, the first `extra` logical row positions get `base+1` rows and the rest get `base`. Compute global starting rows from prefix sums or the equivalent arithmetic formula.

Use a process topology with upper/lower neighbors. Do not confuse Cartesian rank numbers with world rank numbers when reordering is enabled.

### Two-dimensional blocks

Use a `[py, px]` Cartesian process grid, with `py*px = p`. Divide rows and columns independently, including remainders. An explicit grid must fit within the global grid dimensions; reject zero-sized local blocks. Automatic process-grid selection should also respect this requirement.

Derive neighbors with MPI topology utilities. Physical walls belong only to ranks on the global edge. Columns are not contiguous in row-major storage: implement either column packing or an appropriate MPI datatype, document it and keep it consistent across communication-mode comparisons.

## Halo exchange

Upwinding and the five-point diffusion term use only the four nearest neighbours, so a one-cell halo is sufficient.

The blocking baseline can use `MPI_Sendrecv` with consistent partners and tags. It must not rely on eager buffering to avoid deadlock. Sides on a physical wall have no rank neighbour; fill those ghosts by the boundary rule instead.

For overlap:

1. Post receives into the current field's halos.
2. Post sends from its owned edge values or dedicated send buffers.
3. Fill physical-wall ghosts from the current owned values; this is local work.
4. Update points that do not depend on incoming halos.
5. Complete all sends and receives.
6. Update the remaining owned points exactly once.
7. Swap current and next buffers.

Do not modify send buffers, read incomplete receive buffers or reuse communication resources before the matching requests complete. Small local blocks may have no strict interior; handle them without duplicate edge updates.

Nonblocking calls do not guarantee a performance improvement. MPI progress, message size and available independent work determine whether overlap helps.

## OpenMP and MPI thread support

Begin with `MPI_Init_thread(..., MPI_THREAD_FUNNELED, ...)` and verify support. Only the initializing thread performs MPI calls. OpenMP worksharing distributes destination grid points among threads. Each thread computes the face fluxes it needs for its own destination cells; do not write shared face-flux arrays from several threads. Complete the worksharing region before the master thread swaps buffers or starts dependent work.

Support builds without OpenMP by guarding OpenMP-specific calls with the build definition. Do not add global per-step barriers without a demonstrated need; halo completion and local worksharing must already enforce their own dependencies.

## Diagnostics

Total scalar `M = sum(c(i,j)) * dx * dy` over interior cells. With `noflux` walls `M` is conserved up to rounding. In the parallel solver `M` needs an `MPI_Allreduce`; compute it only at snapshot intervals or after the timed loop, and compare serial and parallel values with a tolerance because the summation order differs.

## Output and metrics

Small validation runs may gather the global matrix. Gathering and file writing happen outside solver-loop timing. Global output contains the `ny x nx` interior cells once, never ghost cells or duplicate process halos. Record the wall type and wall values in the run metadata.

Time the whole loop, including necessary communication, wall-ghost filling and any per-step velocity evaluation. Synchronize before the timed MPI region and report the maximum elapsed time over ranks. Keep initialization, validation and output timings separate. Phase timings under overlap are diagnostic and must not be blindly added together.

## Scope

The project uses CPUs, a five-point stencil with first-order upwind advection, constant diffusivity, a prescribed divergence-free velocity field and rectangular geometry. Solving for the velocity (Navier–Stokes or Boussinesq flow), higher-order advection schemes with wider halos, variable materials, complex geometry, implicit methods, CUDA and a web interface are outside the current plan.

## Follow-up from the advection–diffusion change

This design supersedes the heat-only, node-based design. The following files still describe the old model and need matching updates:

- `docs/experiment_plan.md`: new checks for mass conservation, constant field with flow, rotation and uniform-flow Gaussian tests; the sine-decay test now uses cell centres.
- `config/defaults.json`, `config/README.md`: `kappa`, `flow` and its parameters, `boundary_type`, `pigment_drop` parameters, the new time-step rule and `dx = lx/nx`.
- `experiments/*.json`: flow and wall settings; weak scaling now uses `lx = px*local_nx*dx`, and needs fixed speed bounds as the domain grows.
- `src/*/README.md`, `results/raw/metrics_columns.csv`, `scripts/README.md`, `hpc/job_template.sh`, `report/report.md` and the top-level `README.md`.
