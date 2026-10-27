# Planned experiments

These JSON documents are experiment specifications, not executable scripts. The planned runner is described in [scripts/README.md](../scripts/README.md).

All experiment files use schema version 2 and inherit [defaults.json](../config/defaults.json). Each file declares its type through `kind`. Paths in `defaults` are relative to that experiment file; result paths are relative to the repository working directory.

| File | Kind | Planned use |
| --- | --- | --- |
| [correctness.json](correctness.json) | `correctness` | Explicit small runs, each compared with its own matching serial reference. |
| [strong_scaling.json](strong_scaling.json) | `strong_scaling` | Fixed global problem and increasing ranks at one thread per rank. |
| [weak_scaling.json](weak_scaling.json) | `weak_scaling` | Fixed interior cells per rank and fixed threads per rank; explicit 2D process grids. |
| [strategies.json](strategies.json) | `strategies` | Fixed total CPU cores while comparing rank/thread allocation and strategies. |

## Run expansion

- For correctness, expand only the explicit `runs` list. Merge `problem` and `execution` overrides recursively into the defaults. Generate a serial reference with the same resolved problem for each run.
- For strong scaling, run each `resources` entry with the shared `problem`, `decomposition` and `exchange` settings. Measure the pure serial executable separately for speedup.
- For all performance experiments, resource entries override the corresponding fields in `execution`; a `problem` object inside a resource entry overrides the corresponding problem fields, and the shared `execution.backend` selects the parallel executable. The pure serial reference is a separately generated configuration with one rank and one thread.
- For strategies, explicitly take the Cartesian product of `resources`, `decompositions` and `exchanges`.
- For weak scaling, use each explicit `resources` entry. With `[py, px]`, derive `nx = local_nx*px` and `ny = local_ny*py`. Use fixed `dx` and `dy` from `spacing`, then derive `lx = nx*dx` and `ly = ny*dy` (cell-centred grid). Each resource entry also sets `flow_amplitude` proportional to `lx`, so the double-gyre speed bounds stay fixed, and centres the blob in the grown domain.

Weak scaling intentionally grows the physical domain. With fixed spacing, fixed $\kappa$ and fixed speed bounds, the automatic `dt` and physical duration stay fixed while preserving the per-rank grid shape. It is a throughput experiment, not a comparison of the same physical domain at several resolutions.

The rank/thread resources are a starting plan, not confirmed Atlas or Vanda allocations. Adjust them after checking the available node layout and queue limits. Store the resolved configuration whenever a plan is changed.

`warmup_runs` are executed but not included in aggregate results. Keep all measured repetitions. Report medians and a measure of variation; do not select only the fastest repetition. Use fixed steps for the primary scaling experiments.
