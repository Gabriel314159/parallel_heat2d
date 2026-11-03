# Experiments

This markdown file explains our experimental design, including one correctness check and five performance experiments. Each experiment changes one factor and keeps the others fixed. The measurement protocol is in [docs/experiment_plan.md](../docs/experiment_plan.md).

| # | Experiment | Question | What changes | What stays fixed | JSON File |
| --- | --- | --- | --- | --- | --- |
| 0 | Correctness | Does every parallel configuration replicate the serial answer? | ranks, threads, decomposition, exchange | small grids; nothing is timed | [correctness.json](correctness.json) |
| 1 | Serial vs MPI (strong scaling over ranks) | How much faster do more ranks make the same problem? | ranks: 1, 2, 4, …, 256 | $P$, one thread per rank, 2D, blocking | [strong_scaling_mpi.json](strong_scaling_mpi.json) |
| 2 | Serial vs OpenMP (strong scaling over threads) | How much faster do more threads make the same problem? | threads: 1, 2, 4, …, 32, 36, 64, 72 (one Vanda node) | $P$, one rank | [strong_scaling_openmp.json](strong_scaling_openmp.json) |
| 3 | Pure MPI vs MPI + OpenMP | With 256 cores, how should they be split between ranks and threads? | ranks × threads: 256×1, 128×2, 64×4, …, 4×64; then 1 or 2 GPUs | $P$, 256 cores for the CPU layouts, 2D, blocking | [hybrid_layouts.json](hybrid_layouts.json), [hybrid_gpu.json](hybrid_gpu.json) |
| 4 | Communication strategy | Do the decomposition and the exchange mode matter? | rows or 2D; blocking or overlap | $P$, 256 ranks × 1 thread | [communication.json](communication.json) |
| 5 | Weak scaling | Does the time stay the same when the work per core stays the same? | cores = ranks: 1, 4, 16, 64, 256 (the domain grows) | 256 × 256 cells per rank, one thread per rank, 2D, blocking | [weak_scaling.json](weak_scaling.json) |

$P$ is the fixed problem specified in [benchmark_problem.json](benchmark_problem.json): 1024 × 1024 cells on the 1 × 1 basin and 100 steps, with the flow starting at rest.

The course asks for scaling tests up to at least 256 cores: experiments 1 and 5 go that far, experiments 3 and 4 compare their variants on 256 cores, and experiment 2 stops at one node, the limit of OpenMP. Every run binds one thread to each physical core and uses no hardware threads (no hyper-threading); see "Binding" in [hpc/README.md](../hpc/README.md). Every scaling result is discussed for the whole computation first, then for its dominant part, the Poisson solve ([docs/experiment_plan.md](../docs/experiment_plan.md), "The whole computation and its dominant part").

## How the experiments fit together

- **Strong scaling is experiments 1 and 2.** Both hold $P$ fixed and add resources: first ranks (the serial vs MPI comparison, one thread per rank), then threads (the serial vs OpenMP comparison, one rank).
- **Experiments 3 and 4 were one file.** `strategies.json` crossed three factors (layout, decomposition and exchange: 20 configurations). Now experiment 3 changes only the layout and experiment 4 changes only the decomposition and the exchange. Interactions between the factors are not measured; if experiment 3 or 4 shows an effect that may depend on the other factor, add a file that crosses them.
- **The baseline is 2D decomposition with blocking exchange** in experiments 1, 2, 3 and 5. Experiment 4 shows whether another choice is better.
- **Experiment 5 is the only one in which the problem grows.** It starts from one core (one rank, one thread; the serial program runs the same problem) and keeps one thread per rank, so the number of cores is the number of ranks.
- **What experiment 4 can show.** `--exchange` selects how both halo exchanges run: the pigment's, once per step, and the conjugate gradient's, once per iteration. The two global sums of every iteration always wait, so overlap can save at most the time of the CG's halo exchanges, which every run reports (`poisson_halo_s`). See the plan.

## Order and repeated runs

Run experiment 0 first, then 1 and 2 (they give the serial reference and the one-core baselines), then 3, 4 and 5. Experiment 0 runs with `scripts/check_correctness.py`, and every other experiment with `scripts/run_benchmarks.py --config experiments/<name>.json` (`--dry-run` shows its configurations and commands first; see `scripts/README.md`).

Some configurations appear in more than one file: both serial references and the 1-rank, 1-thread run (experiments 1 and 2), and 256 ranks × 1 thread with 2D and blocking (1, 3 and 4). Of the 33 configurations listed in experiments 1 to 4 (CPU), 28 are distinct; the GPU runs of experiment 3 add 2. Each file stays runnable on its own, so a repeated configuration is measured again; if cluster time is short, run the files in one session and use the same rows for the repeats.

## The shared problem

[benchmark_problem.json](benchmark_problem.json) holds $P$ once. Experiments 1 to 4 inherit it, so they solve the identical problem and one serial reference serves all four; experiment 5 inherits its flow and pigment settings and replaces the grid. The values are written out in full, so that a change of `config/defaults.json` cannot change a benchmark.

The grid is fixed at 1024 × 1024 cells; the number of steps is provisional. A serial run of 100 steps takes about 8 minutes on a laptop. Time one run on Vanda, then adjust `steps` in this one file.

## File format

- Every file has `schema_version` 5, the same number as [config/defaults.json](../config/defaults.json). 
- `defaults` names the file to inherit, relative to the file that names it. The experiments inherit `benchmark_problem.json`, which inherits `config/defaults.json`; the chain is resolved from the inside out (`resolve_base` in [scripts/heat2d_tools.py](../scripts/heat2d_tools.py)). Result paths are relative to the repository working directory.
- Each experiment writes its raw measurements to `results/raw/<file name>.csv` and its other files to `outputs/<file name>/`.

## Run expansion

- For correctness (`correctness`), expand only the explicit `runs` list. Merge `problem` and `execution` overrides recursively into the defaults. Generate a serial reference with the same resolved problem for each run.
- For strong scaling (`strong_scaling`), run each `resources` entry with the shared `decomposition` and `exchange`. With `include_serial_reference`, also run both one-core references on the same problem: the pure serial executable `heat_serial` and `heat_serial_omp`, the same solver compiled with the OpenMP pragmas and run on one thread (see `docs/experiment_plan.md`, "Two serial references").
- For the fixed-core comparisons (`strategies`), explicitly take the Cartesian product of `resources`, `decompositions` and `exchanges`. Every resource entry must use `total_cores` cores. Experiment 3 has one decomposition and one exchange, so the product is its seven layouts; experiment 4 has one layout, so the product is its four variants.
- For weak scaling (`weak_scaling`), use each explicit `resources` entry. With `[py, px]`, derive `nx = local_nx*px` and `ny = local_ny*py`. Use fixed `dx` and `dy` from `spacing`, then derive `lx = nx*dx` and `ly = ny*dy` (cell-centred grid). Each resource entry centres the blob in the grown domain. The flow settings and the speed estimate are fixed, so `dt` is the same for every run. With `include_serial_reference`, also run both serial references on the problem of the first entry (one rank, one thread).
- For all performance experiments, resource entries override the corresponding fields in `execution`; a `problem` object inside a resource entry overrides the corresponding problem fields, and the shared `execution.backend` selects the parallel executable. The pure serial reference is a separately generated configuration with one rank and one thread.

Weak scaling intentionally grows the physical domain. With fixed spacing, fixed $\kappa$ and fixed speed bounds, the automatic `dt` and physical duration stay fixed while preserving the per-rank grid shape. It is a throughput experiment, not a comparison of the same physical domain at several resolutions.

The resources follow Vanda's nodes: 72 cores in two sockets of 36 ([hpc/README.md](../hpc/README.md)). Experiment 2 goes up to one node, 72 threads. Runs on several nodes use 64 cores of each node, 32 per socket, so the 256-core runs of experiments 1, 3, 4 and 5 take four nodes; they need the charged queue `batch_cpu`, as the free queue stops at 160 cores. Store the resolved configuration whenever a plan is changed.

`warmup_runs` are executed but not included in aggregate results. Keep all measured repetitions. Report medians and a measure of variation; do not select only the fastest repetition. Use fixed steps for the primary scaling experiments.
