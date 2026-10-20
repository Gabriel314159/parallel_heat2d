# Configuration conventions

These JSON files define the planned interface. No program reads them yet. The future Python runner will load them and invoke the C++ CLI. Individual JSON paths and runtime output paths are resolved as described below.

## Problem parameters

| JSON field | Planned CLI | Meaning |
| --- | --- | --- |
| `problem.nx`, `problem.ny` | `--nx`, `--ny` | Global numbers of interior columns and rows. |
| `problem.lx`, `problem.ly` | `--lx`, `--ly` | Physical plate dimensions. |
| `problem.steps` | `--steps` | Fixed number of explicit time steps. |
| `problem.alpha` | `--alpha` | Positive constant diffusivity. |
| `problem.cfl` | `--cfl` | Fraction of the maximum stable time step; must be in `(0,1]`. |
| `problem.dt` | `--dt` | Optional positive time step. `null` selects automatic time stepping. |
| `problem.case` | `--case` | `heated_top`, `constant`, or `sine_decay`. |
| `problem.initial_temperature` | `--initial-temperature` | Initial interior temperature for heated-top and constant cases. |
| `problem.boundaries.*` | `--top`, `--bottom`, `--left`, `--right` | Fixed physical boundary values. |
| `problem.amplitude` | `--amplitude` | Amplitude for `sine_decay`; defaults to 1.0. |

Explicit `dt` takes precedence over automatic `cfl` time-step selection, but still must satisfy `rx + ry <= 0.5`. Derive `dx` and `dy` from global dimensions. Never derive them from local subdomain sizes.

`constant` requires the initial temperature and all four boundaries to agree. `sine_decay` requires zero values on all four boundaries and initializes `A*sin(pi*x/lx)*sin(pi*y/ly)`; it does not use `initial_temperature`.

The top row of exported matrices corresponds to the physical top boundary (`y=ly`). For the heated-top example, the four corners use the left/right boundary values. Five-point updates do not directly read the corner points, but this convention keeps full-field comparisons deterministic.

## Execution parameters

| JSON field | Planned setting | Meaning |
| --- | --- | --- |
| `execution.backend` | Executable selection | `serial` or `parallel`. |
| `execution.mpi_ranks` | MPI launch option | Number of MPI ranks. |
| `execution.omp_threads` | `OMP_NUM_THREADS` | Threads per rank; also verify the actual number used. |
| `execution.decomposition` | `--decomp` | `rows` or `2d`. |
| `execution.exchange` | `--exchange` | `blocking` or `overlap`. |
| `execution.process_grid` | `--py`, `--px` | Optional `[py, px]`; `py*px` must equal the number of ranks. |

`serial` requires one rank and one thread. MPI must be launched by the appropriate local or HPC launcher. An unspecified 2D process grid will be chosen using MPI topology utilities. Reject grids that would leave any rank with zero interior rows or columns. Handle local dimensions of one or two cells correctly, even when no strictly interior region remains for overlap.

## Output parameters

- `output.directory`: relative to the repository working directory.
- `output.write_field`: whether to write the final temperature matrix.
- `output.snapshot_interval`: steps between snapshots; zero disables intermediate output.
- Planned `--output PATH`: writes a selected final field. Create parent directories as needed.

Validation dumps are small. Large benchmark runs should normally disable snapshots and full-field output. Field matrices are CSV files with `(ny+2)` rows and `(nx+2)` columns, including physical boundaries and excluding process halos. They contain numeric values only, with enough precision for double-precision round trips.

## Defaults and overrides

An experiment's `defaults` path is relative to that experiment JSON file. Load the defaults, then recursively merge object fields from the experiment and the individual run. Scalars and arrays replace previous values. The runner must not pass experiment metadata such as `schema_version`, `checks` or `repetitions` to the solver.

Experiment files define explicit resource lists or explicit run lists. Do not silently interpret lists of ranks and threads as a Cartesian product unless that file says so. Save the fully resolved configuration with each measurement.
