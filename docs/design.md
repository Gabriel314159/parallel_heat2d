# Design notes

Status: implemented CPU design. For completed local checks, see [verification.md](verification.md). HPC measurements remain to be collected.

## Model and indexing

Use constant-coefficient diffusion on a rectangular plate, with fixed Dirichlet boundary conditions and no source term. Start with an explicit forward-Euler/five-point finite-difference method. All implementations solve the same discrete problem.

- `nx`: global interior columns, along x.
- `ny`: global interior rows, along y.
- `dx = lx/(nx+1)`, `dy = ly/(ny+1)`.
- Grid access uses `(row, column)`; top-down row ordering is used in field files.
- For global interior row `i=1..ny`, use `y = ly-i*dy`; for column `j=1..nx`, use `x = j*dx`.
- The row at `i=0` is the top physical edge and the row at `i=ny+1` is the bottom physical edge.

The heated-top example has incompatible side temperatures at the corner points. Assign the four corners from the left/right edge values in every implementation. The stencil for an interior point never directly uses a corner value.

## Time step

Define `rx = alpha*dt/dx^2` and `ry = alpha*dt/dy^2`. Validate positive dimensions, positive `alpha`, positive `dt` and `rx+ry <= 0.5`. For automatic time stepping, `cfl` is a safety fraction of the largest stable time step:

```text
dt = cfl / (2*alpha*(1/dx^2 + 1/dy^2))
```

Use fixed steps for primary scaling tests. Comparisons at different grid resolutions over the same physical duration must recompute stable time steps and the corresponding step counts; record the actual final time.

## Data storage and shared computation

Use contiguous `std::vector<double>` storage. Each local field includes a one-cell halo on every side, and each solver owns a current and next field. The row stride is `local_nx+2`.

The shared stencil computes a single destination value and contains no MPI calls, OpenMP regions, allocation or output. The serial and parallel solvers control the loops. This separation makes formula changes explicit and keeps the serial executable free from MPI/OpenMP dependencies.

Initialize both buffers' physical edges. Initialize interiors by global coordinates so a change in decomposition cannot change the initial condition. Initialize received halos through exchange before their values are used.

## Domain decomposition

### Rows

Distribute `ny` rows over the ranks. With `base=ny/p` and `extra=ny%p`, the first `extra` logical row positions get `base+1` rows and the rest get `base`. Compute global starting rows from prefix sums or the equivalent arithmetic formula.

Use a process topology with upper/lower neighbors. Do not confuse Cartesian rank numbers with world rank numbers when reordering is enabled.

### Two-dimensional blocks

Use a `[py, px]` Cartesian process grid, with `py*px = p`. Divide rows and columns independently, including remainders. An explicit grid must fit within the interior grid dimensions; reject zero-sized local blocks. Automatic process-grid selection should also respect this requirement.

Derive neighbors with MPI topology utilities. Physical boundaries belong only to ranks on the global edge. Columns are not contiguous in row-major storage: implement either column packing or an appropriate MPI datatype, document it and keep it consistent across communication-mode comparisons.

## Halo exchange

The blocking baseline can use `MPI_Sendrecv` with consistent partners and tags. It must not rely on eager buffering to avoid deadlock. The physical exterior has no rank neighbor; preserve its fixed boundary values.

For overlap:

1. Post receives into the current field's halos.
2. Post sends from its owned edge values or dedicated send buffers.
3. Update points that do not depend on incoming halos.
4. Complete all sends and receives.
5. Update the remaining owned points exactly once.
6. Swap current and next buffers.

Do not modify send buffers, read incomplete receive buffers or reuse communication resources before the matching requests complete. Small local blocks may have no strict interior; handle them without duplicate edge updates.

Nonblocking calls do not guarantee a performance improvement. MPI progress, message size and available independent work determine whether overlap helps.

## OpenMP and MPI thread support

Begin with `MPI_Init_thread(..., MPI_THREAD_FUNNELED, ...)` and verify support. Only the initializing thread performs MPI calls. OpenMP worksharing distributes destination grid points among threads. Complete the worksharing region before the master thread swaps buffers or starts dependent work.

Support builds without OpenMP by guarding OpenMP-specific calls with the build definition. Do not add global per-step barriers without a demonstrated need; halo completion and local worksharing must already enforce their own dependencies.

## Output and metrics

Small validation runs may gather the global matrix. Gathering and file writing happen outside solver-loop timing. Global output contains physical boundaries once and never contains duplicate process halos.

Time the whole loop, including necessary communication. Synchronize before the timed MPI region and report the maximum elapsed time over ranks. Keep initialization, validation and output timings separate. Phase timings under overlap are diagnostic and must not be blindly added together.

## Scope

The initial project uses CPUs, fixed boundaries, constant diffusivity and a five-point stencil. CUDA, variable materials, complex geometry, implicit methods and a web interface are outside the current plan.
