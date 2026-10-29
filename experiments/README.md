# Experiments

This markdown file explains our experimental design, including one correctness check and five performance experiments. Each experiment changes one factor and keeps the others fixed. The measurement protocol is in [docs/experiment_plan.md](../docs/experiment_plan.md).

| # | Experiment | Question | What changes | What stays fixed | JSON File |
| --- | --- | --- | --- | --- | --- |
| 0 | Correctness | Does every parallel configuration replicate the serial answer? | ranks, threads, decomposition, exchange | small grids; nothing is timed | [correctness.json](correctness.json) |
| 1 | Serial vs MPI (strong scaling over ranks) | How much faster do more ranks make the same problem? | ranks: 1, 2, 4, 8, 16 | $P$, one thread per rank, 2D, blocking | [strong_scaling_mpi.json](strong_scaling_mpi.json) |
| 2 | Serial vs OpenMP (strong scaling over threads) | How much faster do more threads make the same problem? | threads: 1, 2, 4, 8, 16 | $P$, one rank | [strong_scaling_openmp.json](strong_scaling_openmp.json) |
| 3 | Pure MPI vs MPI + OpenMP | With 16 cores, how should they be split between ranks and threads? | ranks × threads: 16×1, 8×2, 4×4, 2×8, 1×16 | $P$, 16 cores, 2D, blocking | [hybrid_layouts.json](hybrid_layouts.json) |
| 4 | Communication strategy | Do the decomposition and the exchange mode matter? | rows or 2D; blocking or overlap | $P$, 16 ranks × 1 thread | [communication.json](communication.json) |
| 5 | Weak scaling | Does the time stay the same when the work per rank stays the same? | ranks: 1, 4, 16 (the domain grows) | 512 × 512 cells and 4 threads per rank | [weak_scaling.json](weak_scaling.json) |

$P$ is the fixed problem specified in [benchmark_problem.json](benchmark_problem.json): 512 × 1024 cells, 1×2 rectangular boundary,  and 100 steps, with the flow starting at rest.

## How the experiments fit together

- **Strong scaling is experiments 1 and 2.** Both hold $P$ fixed and add resources: first ranks (the serial vs MPI comparison, one thread per rank), then threads (the serial vs OpenMP comparison, one rank).
- **Experiments 3 and 4 were one file.** `strategies.json` crossed three factors (layout, decomposition and exchange: 20 configurations). Now experiment 3 changes only the layout and experiment 4 changes only the decomposition and the exchange. Interactions between the factors are not measured; if experiment 3 or 4 shows an effect that may depend on the other factor, add a file that crosses them.
- **The baseline is 2D decomposition with blocking exchange** in experiments 1, 2, 3 and 5. Experiment 4 shows whether another choice is better.
- **Experiment 5 is the only one in which the problem grows.**
- **What experiment 4 can show.** `--exchange` selects how both halo exchanges run: the pigment's, once per step, and the conjugate gradient's, once per iteration. The two global sums of every iteration always wait, so overlap can save at most the time of the CG's halo exchanges, which every run reports (`poisson_halo_s`). See the plan.

## Order and repeated runs

Run experiment 0 first, then 1 and 2 (they give the serial reference and the one-core baselines), then 3, 4 and 5.

Some configurations appear in more than one file: the serial run and the 1-rank, 1-thread run (experiments 1 and 2), 16 ranks × 1 thread with 2D and blocking (1, 3 and 4), and 1 rank × 16 threads (2 and 3). Of the 21 configurations listed in experiments 1 to 4, 16 are distinct. Each file stays runnable on its own, so a repeated configuration is measured again; if cluster time is short, run the files in one session and use the same rows for the repeats.

## The shared problem

[benchmark_problem.json](benchmark_problem.json) holds $P$ once. Experiments 1 to 4 inherit it, so they solve the identical problem and one serial reference serves all four; experiment 5 inherits its flow and pigment settings and replaces the grid. The values are written out in full, so that a change of `config/defaults.json` cannot change a benchmark.

The grid and the number of steps are provisional. They are chosen so that a serial run takes a few minutes on a laptop (the first step took 1.7 s). Time one run on Atlas or Vanda, then adjust `nx`, `ny` and `steps` in this one file.

## File format

- Every file has `schema_version` 5, the same number as [config/defaults.json](../config/defaults.json). 
- `defaults` names the file to inherit, relative to the file that names it. The experiments inherit `benchmark_problem.json`, which inherits `config/defaults.json`; the chain is resolved from the inside out (`resolve_base` in [scripts/heat2d_tools.py](../scripts/heat2d_tools.py)). Result paths are relative to the repository working directory.
- Each experiment writes its raw measurements to `results/raw/<file name>.csv` and its other files to `outputs/<file name>/`.

## Run expansion

- For correctness (`correctness`), expand only the explicit `runs` list. Merge `problem` and `execution` overrides recursively into the defaults. Generate a serial reference with the same resolved problem for each run.
- For strong scaling (`strong_scaling`), run each `resources` entry with the shared `decomposition` and `exchange`. With `include_serial_reference`, also run the pure serial executable on the same problem, as the reference for the speedup.
- For the fixed-core comparisons (`strategies`), explicitly take the Cartesian product of `resources`, `decompositions` and `exchanges`. Every resource entry must use `total_cores` cores. Experiment 3 has one decomposition and one exchange, so the product is its five layouts; experiment 4 has one layout, so the product is its four variants.
- For weak scaling (`weak_scaling`), use each explicit `resources` entry. With `[py, px]`, derive `nx = local_nx*px` and `ny = local_ny*py`. Use fixed `dx` and `dy` from `spacing`, then derive `lx = nx*dx` and `ly = ny*dy` (cell-centred grid). Each resource entry centres the blob in the grown domain. The flow settings and the speed estimate are fixed, so `dt` is the same for every run.
- For all performance experiments, resource entries override the corresponding fields in `execution`; a `problem` object inside a resource entry overrides the corresponding problem fields, and the shared `execution.backend` selects the parallel executable. The pure serial reference is a separately generated configuration with one rank and one thread.

Weak scaling intentionally grows the physical domain. With fixed spacing, fixed $\kappa$ and fixed speed bounds, the automatic `dt` and physical duration stay fixed while preserving the per-rank grid shape. It is a throughput experiment, not a comparison of the same physical domain at several resolutions.

The rank/thread resources are a starting plan, not confirmed Atlas or Vanda allocations. Adjust them after checking the available node layout and queue limits. Store the resolved configuration whenever a plan is changed.

`warmup_runs` are executed but not included in aggregate results. Keep all measured repetitions. Report medians and a measure of variation; do not select only the fastest repetition. Use fixed steps for the primary scaling experiments.
