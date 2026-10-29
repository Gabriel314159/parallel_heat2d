# Parallel Heat2D

A course project on a parallel two-dimensional solver for pigment carried by a wind-driven flow in a closed basin, using C++, MPI and OpenMP.

The same numerical method runs in serial and in parallel, to compare how domain decomposition, halo communication and thread allocation affect correctness and performance.

## Current status

Implemented and tested:

- Shared components ([src/common/](src/README.md)): command line and validation, grid storage, the flow (vorticity equation with the Arakawa Jacobian and a conjugate gradient Poisson solver), initial conditions, no-flux walls, the upwind flux stencil, and field and snapshot output.
- The serial solver `heat_serial` and the MPI + OpenMP solver `heat_parallel` (row and 2D decompositions, blocking and overlapping halo exchange). On every correctness case, `heat_parallel` agrees with the serial result to within about $5\times10^{-14}$, flow included. The match is not bit for bit, because the global sums of the Poisson solver add the ranks' contributions in a different order.
- Unit tests ([tests/](tests/), run with `ctest`) and [scripts/check_correctness.py](scripts/check_correctness.py): mass conservation, constant fields under the flow, the exact discrete cosine-decay solution, and serial/parallel comparisons.
- [scripts/plot_field.py](scripts/plot_field.py) and [scripts/animate_field.py](scripts/animate_field.py) for pictures and animations of the pigment.

Not implemented yet:

- The benchmark runner and scaling plots (`scripts/run_benchmarks.py`, `scripts/plot_results.py`).
- Measurements on Atlas or Vanda.

## Problem

A pigment with concentration $c(x,y,t)$ is carried by a flow $\mathbf{u}$ and diffuses in a closed basin $[0,l_x]\times[0,l_y]$:

$$
\frac{\partial c}{\partial t} + \nabla\cdot(\mathbf{u}c) = \kappa\nabla^2 c .
$$

The flow is computed at the same time, from the two-dimensional Navier–Stokes equations in vorticity form, driven by the wind and damped by viscosity:

$$
\frac{\partial \zeta}{\partial t} + \mathbf{u}\cdot\nabla\zeta = F(y,t) + \nu\nabla^2\zeta,
\qquad
\nabla^2\psi = -\zeta,
\qquad
u = \frac{\partial\psi}{\partial y},\quad v = -\frac{\partial\psi}{\partial x}.
$$

- The wind $F = F_0\sin(2\pi y/l_y)$ (`--wind steady`) drives two gyres; `--wind periodic` lets its strength oscillate, $F_0(1+\alpha\sin(\omega t))\sin(2\pi y/l_y)$.
- The pigment walls are no-flux walls, so the total pigment is conserved; the flow walls are free-slip ($\psi = \zeta = 0$).
- The pigment does not affect the flow. With no wind, the flow stays at rest and the pigment obeys the heat equation, which gives an analytical reference.

Numerics: a cell-centred finite-volume grid for the pigment (first-order upwind advection, central diffusion, forward Euler) and the cell corners for $\psi$ and $\zeta$ (Arakawa Jacobian, third-order Adams–Bashforth, and a conjugate gradient Poisson solver). The velocities on the cell faces come from $\psi$ and are exactly divergence-free. The time step is fixed for a run, derived from an estimate of the largest speed (`--max-speed`), and the actual speeds are checked during the run.

The full specification is in [docs/design.md](docs/design.md).

## Build and test

Requirements: CMake 3.16 or newer and a C++17 compiler (for example GCC 9+ or Clang 10+). `heat_parallel` also needs an MPI library and, for threads, OpenMP. Without MPI, only `heat_serial` is built, with a warning; without OpenMP, `heat_parallel` runs one thread per rank.

- Ubuntu/Debian: `sudo apt install cmake g++ libopenmpi-dev openmpi-bin`
- macOS: `brew install cmake open-mpi libomp` (CMake finds Homebrew's libomp itself)
- Clusters: load the compiler and MPI modules first, e.g. `module load ...`

```bash
cmake -S . -B build
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The build type defaults to Release (`-O3`). The `correctness_parallel` test launches up to 6 MPI ranks; on a machine with fewer cores, Open MPI needs `-DHEAT2D_MPIEXEC_FLAGS=--oversubscribe`. `-DHEAT2D_NATIVE=ON` also optimises for the build machine's CPU (`-march=native`); such a binary may not run on other machines.

## Run the solvers

```bash
./build/heat_serial --help
./build/heat_serial --nx 256 --ny 256 --steps 10000 --output outputs/final.csv

# 4 MPI ranks x 2 OpenMP threads, 2D process grid, overlapping halo exchange
OMP_NUM_THREADS=2 mpiexec -n 4 ./build/heat_parallel --nx 256 --ny 256 --steps 10000 \
    --decomp 2d --exchange overlap --output outputs/final_parallel.csv
```

Both solvers take the same options; `heat_parallel` adds `--decomp rows|2d`, `--exchange blocking|overlap`, `--py N --px N` for a fixed process grid, and `--threads N` for the OpenMP threads per rank. Each run prints a short report and ends with a `RESULT {...}` line of JSON: all parameters, `dt`, the solver-loop time (the maximum over ranks), the time of the Poisson solves with the parts spent in their halo exchanges and global sums, the relative mass change and the Péclet numbers. With `--output`, the final field is written as CSV (or as binary `.npy` if the path ends in `.npy`) and the same JSON is written next to it. [config/README.md](config/README.md) lists every option.

| Mode | Purpose |
| --- | --- |
| `heat_serial` | Single-threaded reference implementation with no MPI dependency. |
| `heat_parallel --decomp rows` | MPI row decomposition; each rank owns a band of rows. |
| `heat_parallel --decomp 2d` | Two-dimensional MPI process grid (`--py`, `--px`, or chosen automatically). |
| `heat_parallel --exchange blocking` | Halo exchange with `MPI_Sendrecv`, then the update. |
| `heat_parallel --exchange overlap` | Nonblocking halo exchanges: the interior of the pigment update and of the CG's $Ap$ is computed while messages travel. |

Threads per rank: `--threads N`, otherwise `OMP_NUM_THREADS`, otherwise 1 (so a run with neither is MPI-only). `heat_serial` accepts only `--threads 1`. The run report and the `RESULT` line state the ranks, threads and nodes actually used, and on Linux the fewest cores any rank may run on (`cores_per_rank`); the solver warns when a rank has fewer cores than threads, which means the launcher bound it too tightly. If anything fails, every rank prints its message and the job is stopped with `MPI_Abort`; Open MPI adds its own abort banner. A single-rank run of `heat_parallel` still includes MPI setup and is not the pure serial reference.

## Preview the dynamics

The plotting scripts need Python 3.9+ with NumPy and Matplotlib 3.6+ (`python3 -m pip install -r requirements.txt`).

```bash
OMP_NUM_THREADS=2 mpiexec -n 4 ./build/heat_parallel --decomp 2d --exchange overlap \
    --snapshot-interval 500 --snapshot-coarsen 2 --snapshot-dir outputs/preview
python3 scripts/animate_field.py outputs/preview --flow streamlines -o outputs/preview.gif
python3 scripts/animate_field.py outputs/preview --field vorticity -o outputs/vorticity.gif
python3 scripts/plot_field.py outputs/preview/field_025000.npy --flow streamlines --gamma 0.5
```

With the defaults ($256 \times 512$ cells, $\nu = 10^{-3}$, $F_0 = 0.1$), the flow first spins up alone for 100 time units into two gyres with speeds of about 1; then the pigment is released on the line between them and followed for 50,000 steps (about 20 time units), with 101 snapshots. That is about 306,000 steps, each needing about 350 CG iterations for $\psi$. On an 8-core laptop with 4 ranks × 2 threads, a step takes about 0.08 s, so the whole run takes roughly 6–7 hours (estimated from the first 700 steps). Run it on a bigger machine, or use a coarser grid for a quick look (with `--nx 128 --ny 256`, each step is about 8× cheaper). Write `.mp4` instead of `.gif` for smoother colours (needs ffmpeg).

Every snapshot holds the pigment, $\psi$ and the vorticity, so `--field vorticity` and `--field psi` show the flow itself and `--flow streamlines` draws the streamlines of the same moment.

## Larger runs on a bigger machine

Each step costs one pigment update, one vorticity update and one Poisson solve. The CG Poisson solve takes over 99% of the time (measured on $128 \times 256$ cells). The time step is set by the speed estimate, $\Delta t = 0.8 \times 0.5/(U_{\max}(1/\Delta x + 1/\Delta y))$, so it halves when the grid is refined. The CG iterations per step also double (about 180 on $128 \times 256$, 350 on $256 \times 512$), so the cost per time unit grows about 16× per doubling of the grid.

The physics also changes with the parameters. With the steady wind and these viscosities, the flow spins up over a viscous time of about $1/(\nu\,\pi^2(1/l_x^2 + 4/l_y^2))$ (25 at $\nu = 2\times10^{-3}$, 100 at $5\times10^{-4}$) into two steady gyres; in test runs up to $\mathrm{Re}$ of about 4000 the flow was still smooth after 100 time units. Use `--wind periodic` (for example `--wind-alpha 1.5 --wind-period 20`, which reverses the wind) for a strongly time-dependent flow.

[hpc/big_run.sh](hpc/big_run.sh) runs one large case with snapshots, serial or parallel; its settings come from environment variables:

```bash
mkdir -p outputs
# 512 x 1024 cells with the default flow, spin-up 100, then about 40 time units with the pigment
HEAT2D_REPO=$PWD RANKS=16 nohup bash hpc/big_run.sh > outputs/run512.log 2>&1 &
```

- `--snapshot-coarsen K` (`COARSEN`) stores the pigment as K × K block averages (which keep the total mass) and $\psi$ and the vorticity at every K-th corner; K must divide `nx` and `ny`.
- The final pigment field goes to `final.npy`; binary output keeps large grids small.
- `nohup` (or `tmux`) keeps the run going after logging out; the report is at the end of the log file.
- If the flow outgrows `--max-speed` (`MAXSPEED`), the run stops with a message naming the limit; rerun with a larger estimate.
- Animate on the big machine, or copy the snapshot folder back, for example `rsync -av machine:parallel_heat2d/outputs/run512 outputs/`.

## Repository layout

| Path | Contents |
| --- | --- |
| [src/](src/README.md) | C++ sources: shared components (`common/`), serial solver (`serial/`), MPI + OpenMP solver (`parallel/`). |
| [tests/](tests/) | C++ unit tests, run by `ctest`. |
| [scripts/](scripts/README.md) | Correctness checks, plotting and animation; planned benchmark tools. |
| [config/defaults.json](config/defaults.json) | Default numerical and execution settings. |
| [experiments/](experiments/README.md) | The correctness check and the five performance experiments, with the shared benchmark problem. |
| [hpc/](hpc/README.md) | HPC setup notes and job templates. |
| [results/raw/](results/raw/README.md) | Measured benchmark CSV files, when available. |
| [results/figures/](results/figures/README.md) | Selected figures generated from measured data. |
| [docs/design.md](docs/design.md) | Numerical and software design decisions. |
| [docs/experiment_plan.md](docs/experiment_plan.md) | Correctness checks and fair timing conventions. |
| [docs/ai_usage.md](docs/ai_usage.md) | Record of AI assistance and subsequent verification. |
| [docs/references.md](docs/references.md) | Technical references. |
| [report/report.md](report/report.md) | Report outline; no experimental findings are claimed. |

Solver output goes to `outputs/`, which Git ignores.

## Validation and experiments

- Compare full scalar fields from serial and parallel runs of the same problem.
- Check mass conservation, constant fields under the flow, non-divisible partitions and small local blocks.
- Use a decaying cosine mode (flow off) for an independent analytical check and mesh refinement study.
- Run the five performance experiments of [experiments/](experiments/README.md): serial vs MPI (strong scaling), serial vs OpenMP, pure MPI vs MPI + OpenMP, communication strategies and weak scaling.
- Exclude file I/O and validation from solver-loop timing, and report them separately where useful.
- Use the maximum elapsed time across MPI ranks as parallel runtime.
- Preserve raw timings, environment details and the exact commit used.

See [experiments/README.md](experiments/README.md) for configuration semantics and [docs/experiment_plan.md](docs/experiment_plan.md) for the measurement protocol.

## Development order

1. Implement shared parameters, grids, initialization and the update formula. (Done.)
2. Implement and validate `heat_serial`. (Done.)
3. Implement row decomposition and blocking MPI halo exchange. (Done.)
4. Add OpenMP computation and verify hybrid configurations. (Done.)
5. Add two-dimensional decomposition, including column communication. (Done.)
6. Add nonblocking exchange and computation overlap. (Done.)
7. Implement the experiment utilities, run local checks and collect HPC measurements.
8. Write conclusions from the measured data and finish the report.

## Submission scope

The assignment asks for a PDF report, source code with running instructions, measured runtime figures from Atlas or Vanda, correctness validation and an explanation of AI use. The supplied assignment specifies submission before 13 November 2026 and mentions an optional progress report before 30 October. Confirm the precise cutoff time on Canvas.
