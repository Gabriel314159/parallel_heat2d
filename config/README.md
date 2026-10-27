# Configuration conventions

These JSON files define the planned interface. No program reads them yet. The future Python runner will load them and invoke the C++ CLI. Individual JSON paths and runtime output paths are resolved as described below.

## Problem parameters

| JSON field | Planned CLI | Meaning |
| --- | --- | --- |
| `problem.nx`, `problem.ny` | `--nx`, `--ny` | Global numbers of cells along x (columns) and y (rows). |
| `problem.lx`, `problem.ly` | `--lx`, `--ly` | Physical domain dimensions. |
| `problem.steps` | `--steps` | Fixed number of explicit time steps. |
| `problem.kappa` | `--kappa` | Positive constant diffusivity $\kappa$. |
| `problem.cfl` | `--cfl` | Safety fraction of the maximum stable time step; must be in `(0,1]`. |
| `problem.dt` | `--dt` | Optional positive time step. `null` selects automatic time stepping. |
| `problem.flow_amplitude` | `--flow-amplitude` | Double-gyre amplitude $A$. `0` switches the flow off (heat equation). |
| `problem.flow_epsilon` | `--flow-epsilon` | Double-gyre oscillation strength $\varepsilon \ge 0$. `0` gives the steady double gyre. |
| `problem.flow_period` | `--flow-period` | Double-gyre period $T > 0$; $\omega = 2\pi/T$. |
| `problem.case` | `--case` | `pigment_drop`, `constant` or `cosine_decay`. |
| `problem.background` | `--background` | Background value $c_0$ for `pigment_drop` and `constant`. |
| `problem.amplitude` | `--amplitude` | Initial-field amplitude $A_0$ for `pigment_drop` and `cosine_decay`. |
| `problem.blob_centre` | `--blob-x`, `--blob-y` | Blob centre $[x_0, y_0]$ for `pigment_drop`, in physical coordinates. |
| `problem.blob_width` | `--blob-width` | Gaussian width $\sigma > 0$ for `pigment_drop`. |

All four walls are no-flux walls; there is no boundary-type option and no wall values.

The grid is cell-centred: $\Delta x = l_x/n_x$ and $\Delta y = l_y/n_y$. Derive them from global dimensions, never from local subdomain sizes.

Time step: with $u_{\max} = \lvert A\rvert\pi/l_y$, $v_{\max} = 2\lvert A\rvert\pi(1+2\varepsilon)/l_x$,

$$
a = \frac{u_{\max}}{\Delta x} + \frac{v_{\max}}{\Delta y},
\qquad
d = 2\kappa\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right).
$$

Automatic stepping uses $\Delta t = \mathrm{cfl}/(a+d)$. An explicit `dt` takes precedence over `cfl` but must still satisfy $\Delta t\,(a+d) \le 1$. See `docs/design.md` for the derivation.

Case requirements:

- `pigment_drop`: initializes $c_0 + A_0\exp\!\big(-((x-x_0)^2+(y-y_0)^2)/(2\sigma^2)\big)$. The default advection–diffusion case; the total mass is conserved.
- `constant`: initializes $c_0$. The field must stay constant under the flow.
- `cosine_decay`: initializes $A_0\cos(\pi x/l_x)\cos(\pi y/l_y)$. Requires `flow_amplitude = 0`. Analytical reference; see `docs/design.md`.

The top row of exported matrices corresponds to the cell row next to the top wall (`y=ly`).

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
- `output.write_field`: whether to write the final scalar field.
- `output.snapshot_interval`: steps between snapshots; zero disables intermediate output.
- Planned `--output PATH`: writes a selected final field. Create parent directories as needed.

Validation dumps are small. Large benchmark runs should normally disable snapshots and full-field output. Field matrices are CSV files with `ny` rows and `nx` columns of interior cell values, excluding ghost cells and process halos. They contain numeric values only, with enough precision for double-precision round trips.

## Defaults and overrides

An experiment's `defaults` path is relative to that experiment JSON file. Load the defaults, then recursively merge object fields from the experiment and the individual run. Scalars and arrays replace previous values. The runner must not pass experiment metadata such as `schema_version`, `checks` or `repetitions` to the solver.

Experiment files define explicit resource lists or explicit run lists. Do not silently interpret lists of ranks and threads as a Cartesian product unless that file says so. Save the fully resolved configuration with each measurement.
