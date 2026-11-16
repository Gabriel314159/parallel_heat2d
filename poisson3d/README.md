# Standalone 3D Poisson benchmark

This directory holds the 3D Poisson solver used for the 2D/3D Poisson comparison of the report. It is a benchmark of its own, separate from the flow and pigment models: the 3D pigment model (model A, `src/common3d`, `src/serial3d`, `src/parallel3d`) still solves the 2D Poisson equation of the flow.

The program solves $-\Delta\phi=f$ on a box $[0,l_x]\times[0,l_y]\times[0,l_z]$ with homogeneous Dirichlet conditions $\phi=0$ on all six faces. The grid parameters $(n_x,n_y,n_z)$ are the numbers of interior unknowns (vertices), with spacing $h_x=l_x/(n_x+1)$ and so on. The seven-point operator is applied matrix-free, and the system is solved with unpreconditioned conjugate gradients (CG), starting from zero. MPI splits the grid in $x$ and $y$ (a $p_x\times p_y$ Cartesian grid, no rank reordering); every rank keeps the whole $z$ direction, as in model A's pencils. Each CG iteration exchanges the halo of the search direction once (blocking, or overlapped with the interior of the matrix–vector product) and does two global sums. OpenMP threads share the loops of a rank, with parallel first-touch allocation. A GPU build runs the CG loops with OpenMP target offload, one GPU per rank; with two ranks the halo goes through host buffers, and the global sums are MPI reductions on the host.

| File | Role |
| --- | --- |
| `poisson3d.cpp` | The solver: options, process layout and binding check, grid and halo datatypes, the seven-point operator, CG, the true-residual check and the JSON result line. Built as `poisson3d_serial` (no OpenMP), `poisson3d_serial_omp`, `poisson3d_mpi` (`-DUSE_MPI`) and, with NVHPC, `poisson3d_omp` and `poisson3d_gpu` (`-DUSE_GPU`). |
| `poisson3d_gpu.hpp` | The CG on the GPU (OpenMP target), with the host-staged halo exchange for two ranks. |
| `build.sh`, `build_gpu.sh` | The builds used for the Vanda measurements: GCC 14.3 + Open MPI 5.0.8 (`-O3 -ffp-contract=off`), and NVHPC 26.1 + HPC-X for the CPU/GPU comparison (`-Mnofma`). |
| `run_experiments.py` | `verify` (validation), `verify-gpu`, and the experiment suites `one-node` and `gpu`; one warm-up and five measured repetitions per configuration, raw CSV plus a JSON record of commands, modules, nodes and binary hashes. |
| `mpi_checked/` | The runner used for the four-node suite: it checks the four-node launch and the binding first, and balances 36 and 144 ranks over the sockets (18 and 24 per socket) as in the 2D weak scaling. See its README. |
| `summarize.py` | Medians and ranges of the measured repetitions, speedups and efficiencies, into one summary CSV. |
| `one_node.pbs`, `four_nodes_checked.pbs`, `two_gpus.pbs`, `one_gpu.pbs` | The Vanda jobs. |

CMake also builds `poisson3d_serial`, `poisson3d_serial_omp` and `poisson3d_mpi` with the rest of the project, and `ctest` runs the serial validation (`poisson3d_verify_serial`).

## Configurations, matched to the 2D experiments

| Experiment | 3D configuration |
| --- | --- |
| Correctness | Refinement against a known continuous solution (second order), an exact discrete solution, point-by-point serial/MPI comparison on uneven blocks with blocking and overlapped exchange, and the one- and two-GPU checks. |
| Serial size | One core, unit cube, 32, 48, 64, 96 and 128 interior points per axis. |
| MPI strong scaling | $128^3$ interior points, 1 to 256 ranks with one thread each; pure serial and the OpenMP-compiled serial program on one thread as references. |
| OpenMP strong scaling | The same $128^3$ grid, one rank with 1, 2, 4, 8, 16, 32, 36, 64 and 72 threads. |
| Hybrid layouts on 256 cores | The same grid, $P\times T$ = $256\times1$, $128\times2$, $64\times4$, $32\times8$, $16\times16$, $8\times32$, $4\times64$. |
| MPI weak scaling | 1, 4, 9, 16, 36, 64, 144 and 256 ranks in one job; $32^3$ interior points per rank, one thread per rank, fixed spacing $1/33$; the box grows in $x$ and $y$ and keeps its depth. The one-rank MPI run is the reference. |
| Decomposition and communication | 256 ranks with one thread each: rows ($1\times256$) or 2D blocks ($16\times16$), blocking or overlapped exchange, on a $256\times1024\times32$ grid, so that rows have 4 rows per rank and an interior to overlap with. |
| CPU/GPU against size | $64^3$, $128^3$, $192^3$, $256^3$; $1\times36$ CPU threads against one A40 with the same layout, and $2\times36$ against two A40. |

Every timed run solves one system with the same right-hand side construction, a zero initial guess and a CG tolerance of $10^{-8}$ on the relative residual. The benchmark right-hand side is a smooth sine plus a deterministic perturbation generated from the global indices, so it does not depend on the decomposition. **The 2D measurements time a whole time-stepping loop (warm-started solves, tolerance $10^{-10}$); the 3D measurements time one cold-started solve.** Compare the two through the time per CG iteration and the shape of the scaling curves, not through iteration counts or total times.

The CPU arrays are first-touched in parallel before timing. The optional controls (the 64 cores of one node as ranks × threads, and 64 ranks on 1, 2 and 4 nodes) are enabled with `qsub -v POISSON3D_CONTROLS=1 ...` or listed with `--include-controls --dry-run`; the measurements in `results/raw/poisson3d/` do not include them.

## Use

From this directory, with C++17 and Python 3.9 or newer. Serial only:

```bash
POISSON3D_SERIAL_ONLY=1 bash build.sh build
python3 run_experiments.py verify --serial-only
./build/poisson3d_serial --nx 32 --ny 32 --nz 32 --rhs benchmark
```

With GCC, Open MPI and OpenMP, `bash build.sh build` builds all CPU programs and `python3 run_experiments.py verify` also compares serial and MPI fields point by point. The verification launches `mpirun` with Open MPI options; `POISSON3D_MPIRUN` replaces the `mpirun` command, but the options stay those of Open MPI.

On Vanda, submit from this directory (the jobs build on the compute node and write to `results/<PBS job ID>/`):

```bash
python3 run_experiments.py one-node --dry-run
qsub -P <project> one_node.pbs
qsub -P <project> four_nodes_checked.pbs
qsub -P <project> two_gpus.pbs
```

`two_gpus.pbs` measures one and two GPUs in one job; `one_gpu.pbs` needs only one GPU and runs the CPU/GPU subset. CPU jobs on four nodes use 64 active physical cores per node; the OpenMP scan of one rank may use all 72 cores of a node. Every run checks the observed ranks, threads, nodes and binding before its row is recorded, and every GPU run checks that the solve ran on the device.

Summarize one set of results with:

```bash
python3 summarize.py one-node.csv four-nodes.csv gpu.csv > summary.csv
```

## Measured results

`results/raw/poisson3d/` holds the Vanda measurements: `one-node.csv` (job 1559053), `four-nodes.csv` (job 1563439, `four_nodes_checked.pbs`) and `gpu.csv` (job 1559087, two A40), each with a JSON record, the four-node launch checks (`preflight-*.json`), the PBS completion records (`jobs.json`) and `summary.csv`. Node labels, GPU UUIDs and home paths were generalized before sharing.

## Columns

- `solve_s`: the initial residual norm, the CG iterations and their communication; for MPI, the maximum over the ranks. Building the right-hand side and the final residual check are outside the timed part.
- `solve_total_s`: `solve_s` plus copying the arrays to and from the GPU; equal to `solve_s` on the CPU. GPU speedups use `solve_total_s`, the time per iteration uses `solve_s` (`solve_s_per_iteration`).
- `halo_s`: time in the halo exchanges, including host staging with two GPUs; with overlap, only posting and waiting. `halo_device_copy_s`: the host/GPU transfers of the two-GPU halo. `reduction_s` and `reduction_min_s`: the largest and smallest time in the global sums over the ranks (the largest includes waiting). These extremes come from different ranks and do not add up.
- `N` $=n_xn_yn_z$; the operator has $7N-2(n_xn_y+n_xn_z+n_yn_z)$ nonzero entries.

Strong-scaling speedups are relative to the pure serial program; weak-scaling efficiency is relative to the one-rank MPI run of the same series. The two GPU ratios are CPU against GPU for the same ranks × threads layout, and one GPU against two (the whole $1\times36$ to $2\times36$ configuration).
