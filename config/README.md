# Configuration conventions

The solvers take command-line flags: the CLI column below, also listed by `heat_serial --help`. The JSON files describe experiments: [scripts/check_correctness.py](../scripts/check_correctness.py) reads them now, and the planned benchmark runner will too. They merge the defaults with the experiment and run overrides, then pass the problem fields to the solver as flags. Individual JSON paths and runtime output paths are resolved as described below.

## Problem parameters

| JSON field | CLI | Meaning |
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

These apply to the parallel solver, `heat_parallel`; `heat_serial` ignores them with a note.

| JSON field | Setting | Meaning |
| --- | --- | --- |
| `execution.backend` | Executable selection | `serial` or `parallel`. |
| `execution.mpi_ranks` | MPI launch option | Number of MPI ranks. |
| `execution.omp_threads` | `OMP_NUM_THREADS` | Threads per rank; also verify the actual number used. |
| `execution.decomposition` | `--decomp` | `rows` or `2d`. |
| `execution.exchange` | `--exchange` | `blocking` or `overlap`. |
| `execution.process_grid` | `--py`, `--px` | Optional `[py, px]`; `py*px` must equal the number of ranks. |

`serial` requires one rank and one thread. MPI must be launched by the appropriate local or HPC launcher. Without `--py`/`--px`, the 2D process grid is chosen with `MPI_Dims_create`, as square as possible, with more process columns than rows when `nx > ny`. Grids that would leave a rank with zero rows or columns are rejected. Local blocks of one or two cells are handled correctly, even when no strictly interior region remains for overlap.

## Output parameters

| JSON field | CLI | Meaning |
| --- | --- | --- |
| `output.directory` | (runner) | Folder for the run's files, relative to the repository working directory. |
| `output.write_field` | `--output PATH` | Write the final field to PATH: `.csv` as text, `.npy` as binary. The run summary is written next to it with the extension `.json`. |
| `output.snapshot_interval` | `--snapshot-interval N` | Write a snapshot every N steps; 0 disables snapshots. |
| | `--snapshot-dir DIR` | Snapshot folder (default `outputs/snapshots`). |
| | `--snapshot-coarsen K` | Store snapshots as K × K block averages (default 1, full size); K must divide `nx` and `ny`. |

Parent folders are created as needed. Field matrices have `ny` rows and `nx` columns of interior cell values, top row first, excluding ghost cells and process halos. CSV files contain numeric values only, with 17 significant digits, enough for double-precision round trips.

A snapshot folder contains `meta.json` (the configuration and `dt`), `index.csv` (one line per snapshot: `step,time,mass,file`) and the fields as `field_NNNNNN.npy`, which `numpy.load` reads. The mass in `index.csv` is that of the full field, also when snapshots are coarsened. Snapshot writing is excluded from the solver-loop timing.

Validation dumps are small. Large benchmark runs should normally disable snapshots and full-field output; for large previews, coarsen the snapshots.

## Defaults and overrides

An experiment's `defaults` path is relative to that experiment JSON file. Load the defaults, then recursively merge object fields from the experiment and the individual run. Scalars and arrays replace previous values. The runner must not pass experiment metadata such as `schema_version`, `checks` or `repetitions` to the solver.

Experiment files define explicit resource lists or explicit run lists. Do not silently interpret lists of ranks and threads as a Cartesian product unless that file says so. Save the fully resolved configuration with each measurement.
