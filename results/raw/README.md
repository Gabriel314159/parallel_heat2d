# Raw measurements

Measured on Vanda (nodes of two Xeon Platinum 8452Y, 72 cores) on 3 and 4 November 2026:

- `strong_scaling_mpi_1node.csv`: experiment 1, both serial references and 1 to 64 ranks on one node.
- `strong_scaling_mpi_4nodes.csv`: experiment 1 at 128 and 256 ranks (64 cores of each node).
- `strong_scaling_openmp.csv`: experiment 2, the references and 1 to 72 threads. Its 9 error rows are a first attempt that the launcher refused (a rank bound across both sockets); they are kept.
- `hybrid_gpu.csv`: experiment 3 with the Poisson solve on one and on two A40 GPUs.
- `communication.csv`: experiment 4, rows or 2D and blocking or overlap, 256 ranks.
- `weak_scaling_1node.csv`: experiment 5, the references and 1, 4, 16 and 64 ranks.
- `weak_scaling_4nodes.csv`: experiment 5 at 256 ranks.
- `correctness_report.json`: experiment 0, 12 cases, no failed check.
- `sessions/<experiment>/<session_id>.json`: the resolved plan and environment of each runner call.

The follow-up runs of 8 November 2026 (`experiments/README.md`, "Follow-up runs"), all at commit b43323f, after the first-touch change and with the column `poisson_sum_min_s`:

- `strong_scaling_openmp_firsttouch.csv`: experiment 2 again, 2 to 72 threads. It ran on a GPU node (`GN-…`, the same CPUs).
- `hybrid_layouts_firsttouch.csv`: experiment 3, the seven layouts of 256 cores from 256 × 1 to 4 × 64, 10 repetitions each.
- `hybrid_node.csv`: the 64 cores of one node as 64 × 1 to 2 × 32 ranks × threads.
- `node_spread.csv`: 64 ranks on 1, 2 and 4 nodes, for the problem of experiment 1 (1024²) and the block per rank of experiment 5 (2048²).
- `weak_scaling_extra_1node.csv` and `weak_scaling_extra_3nodes.csv`: experiment 5 at 9 and 36 ranks, and at 144 ranks.
- `placement_control.csv` and `placement_control_mhz.log`: the control runs on one node (16 and 32 ranks on socket 0 or split, 64 ranks twice) and the CPU clock of every core every 5 s during the second 64-rank run.
- `gpu_sizes_cpu.csv` and `gpu_sizes_gpu.csv`: one rank and two ranks of 36 threads, on the CPU and with one and two A40, at 512², 1024², 2048² and 4096² (100, 100, 20 and 10 steps).

`git_commit` is `unknown` in the rows measured on CPU nodes (`CN-…`), where `git` was not found. Judging by the session times, the four-node runs and the experiment 2 and 5 runs ran at these commits: the four-node runs at 9bafe9e, experiment 5 (one node) at 445dd60, and in experiment 2 the references and error rows at 445dd60 and the threaded rows at 2340aab.

[metrics_columns.csv](metrics_columns.csv) is a header-only format template. It must not be treated as experimental data.

Commit compact raw CSV files and associated resolved configurations/environment records when experiments are run. Each experiment of `experiments/` has its own file, named after its experiment file (for example `strong_scaling_mpi.csv`); `scripts/run_benchmarks.py` appends to it. Experiments 1 and 5 run in two jobs that may run at the same time, so their rows go to `<name>_1node.csv` and `<name>_4nodes.csv` (`--metrics-csv`); join them before the analysis. Keep every measured repetition, its success/failure status and the exact Git commit.

`session_id` names one invocation of the runner; its resolved plan and environment are in `outputs/<experiment>/<session_id>/session.json` (copies of the finished ones in `sessions/`). `phase` is `warmup` or `measured`: warmup rows are kept but left out of every aggregate. `status` is `ok` or `error`, and `error` holds the message of a failed run. `git_commit` ends in `-dirty` when tracked files differed from that commit, and in `-unchecked` when the node has no `git` and the runner only read `.git` (then the files were not compared). A row's `elapsed_s` is the maximum solver-loop time across ranks for a parallel run, or the serial solver-loop time for the serial reference.

`total_cores = mpi_ranks * omp_threads`; In rows with status `ok`, `nodes` and the observed thread count (`omp_threads`) come from the solver's `RESULT` line; error rows hold the requested values; record the placement in session metadata. `cores_per_rank` is the fewest physical cores any rank's threads may run on (0 if the system does not report it); a row with `cores_per_rank` below `omp_threads` had threads sharing cores. `binding` says how the threads were bound to physical cores: `cores` (every thread on a core of its own), `shared` (two threads on one core, for example its two hardware threads), `unpinned` (a thread free to run on several cores) or `unknown`. Only rows with `binding` `cores` are valid timings. `poisson_device` is where the CG loops ran (`cpu`; in the GPU build `gpu`, or `host (no GPU found)`), and `poisson_copy_s` the time of the copies between host and GPU (0 on the CPU). `poisson_s` is the time of the Poisson solves in the timed loop; `poisson_halo_s` and `poisson_sum_s` are the parts of it spent inside the CG's halo exchanges (under overlap, only posting and waiting) and global sums. All three are maxima over the ranks, so they do not add up. `poisson_sum_min_s` is the smallest global-sum time over the ranks, close to the reduction alone (an upper bound on it); the files written before 5 November 2026 do not have this column, so join them with the others by column name. `max_abs_error` is optional for performance-only runs and should be blank when not checked. `mass_error` is the change of the total scalar relative to its initial absolute amount, $\lvert M_{\text{end}}-M_0\rvert / (\Delta x\,\Delta y\sum\lvert c^0_{i,j}\rvert)$, computed outside the timed loop. For non-negative fields this is $\lvert M_{\text{end}}-M_0\rvert/M_0$; unlike that ratio, it is also defined for `cosine_decay`, where $M_0 = 0$. `peclet` and `peclet_grid` are $\mathrm{Pe}$ and $\mathrm{Pe}_{\text{grid}}$ as defined in `docs/design.md`. Analytical discretization errors need a separate validation record.

Do not commit full-resolution field dumps, build directories or temporary log bundles here. Large temporary files belong in ignored output directories or the appropriate HPC filesystem.

The 3D pigment experiments (`experiments/*_3d.json`) write the columns of `metrics_columns_3d.csv`: those of the 2D files without the GPU columns, plus `nz`, `lz`, `flow` (`live` or `frozen`), `pigment_s` (the pigment updates of the timed loop, halo exchanges included), `pigment_halo_s` (those halo exchanges; under overlap, only posting and waiting) and `cell_updates_per_s` ($n_x n_y n_z$ cells per step). Both pigment times are maxima over the ranks.
