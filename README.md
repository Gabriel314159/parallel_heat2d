# Parallel Heat2D

A course project on a parallel two-dimensional advection–diffusion solver (pigment in a double-gyre flow) using C++, MPI and OpenMP.

We will implement the same numerical method in serial and in parallel, then compare how domain decomposition, halo communication and thread allocation affect correctness and performance.

## Current status

Implemented and tested:

- Shared components ([src/common/](src/README.md)): command line and validation, grid storage, the double-gyre velocity, initial conditions, no-flux walls, the upwind flux stencil, and field and snapshot output.
- The serial solver `heat_serial` and the MPI + OpenMP solver `heat_parallel` (row and 2D decompositions, blocking and overlapping halo exchange). On every correctness case, `heat_parallel` reproduces the serial field bit for bit.
- Unit tests ([tests/](tests/), run with `ctest`) and [scripts/check_correctness.py](scripts/check_correctness.py): mass conservation, constant fields under the flow, the exact discrete cosine-decay solution, and serial/parallel comparisons.
- [scripts/plot_field.py](scripts/plot_field.py) and [scripts/animate_field.py](scripts/animate_field.py) for pictures and animations of the pigment.

Not implemented yet:

- The benchmark runner and scaling plots (`scripts/run_benchmarks.py`, `scripts/plot_results.py`).
- Measurements on Atlas or Vanda.

Checks of the original repository scaffold are recorded in [docs/scaffold_checks.md](docs/scaffold_checks.md).

## Problem

We model a passive scalar $c(x,y,t)$, for example pigment concentration, carried by a prescribed flow and diffusing in a closed rectangular box $[0,l_x]\times[0,l_y]$:

$$
\frac{\partial c}{\partial t} + \nabla\cdot(\mathbf{u}c) = \kappa\left(\frac{\partial^2 c}{\partial x^2} + \frac{\partial^2 c}{\partial y^2}\right)
$$

- $\kappa > 0$ is a constant diffusivity.
- $\mathbf{u} = (u, v)$ is the time-periodic **double gyre**, defined by a streamfunction $\psi$ with $u = \partial\psi/\partial y$, $v = -\partial\psi/\partial x$. It is divergence-free and has zero normal velocity on every wall.
- All four walls are **no-flux** walls, so the total amount of scalar is conserved.
- The default case drops a Gaussian blob of pigment into the box and lets the gyres stretch and fold it.
- With flow amplitude $A = 0$ the model reduces to the heat equation with insulated walls, which gives an analytical reference.

The numerical method is a cell-centred finite-volume scheme: first-order upwind advection, central-difference diffusion and forward Euler in time, on a grid with $\Delta x = l_x/n_x$ and $\Delta y = l_y/n_y$. `nx` and `ny` count cells. The time step must satisfy

$$
\Delta t\,(a + d) \le 1,
\qquad
a = \frac{u_{\max}}{\Delta x} + \frac{v_{\max}}{\Delta y},
\qquad
d = 2\kappa\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right),
$$

and automatic stepping uses $\Delta t = \mathrm{cfl}/(a+d)$ with $0 < \mathrm{cfl} \le 1$ (default 0.8).

The full specification, including the face velocities, flux formulas and parallel design, is in [docs/design.md](docs/design.md).

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

Both solvers take the same options; `heat_parallel` adds `--decomp rows|2d`, `--exchange blocking|overlap` and `--py N --px N` for a fixed process grid. Each run prints a short report and ends with a `RESULT {...}` line of JSON: all parameters, `dt`, the solver-loop time (the maximum over ranks), the relative mass change and the Péclet numbers. With `--output`, the final field is written as CSV (or as binary `.npy` if the path ends in `.npy`) and the same JSON is written next to it. [config/README.md](config/README.md) lists every option.

| Mode | Purpose |
| --- | --- |
| `heat_serial` | Single-threaded reference implementation with no MPI dependency. |
| `heat_parallel --decomp rows` | MPI row decomposition; each rank owns a band of rows. |
| `heat_parallel --decomp 2d` | Two-dimensional MPI process grid (`--py`, `--px`, or chosen automatically). |
| `heat_parallel --exchange blocking` | Halo exchange with `MPI_Sendrecv`, then the update. |
| `heat_parallel --exchange overlap` | Nonblocking halo exchange; the interior is updated while messages travel. |

Set `OMP_NUM_THREADS=1` for an MPI-only baseline. A single-rank run of `heat_parallel` still includes MPI setup and is not the pure serial reference.

## Preview the dynamics

The plotting scripts need Python 3.9+ with NumPy and Matplotlib 3.6+ (`python3 -m pip install -r requirements.txt`).

```bash
./build/heat_serial --nx 512 --ny 512 --kappa 0.001 --steps 21200 \
    --snapshot-interval 212 --snapshot-dir outputs/preview
python3 scripts/animate_field.py outputs/preview --velocity streamlines -o outputs/preview.gif
python3 scripts/plot_field.py outputs/preview/field_004240.npy --gamma 0.5 --velocity streamlines
```

This simulates one flow period ($t \approx 10$) on a 512 × 512 grid in well under a minute on one laptop core, and writes 101 snapshots of 2 MB each. The blob starts on the line between the two gyres, is carried up along it, and is then stretched around both gyres; the gyre centres stay pigment-poor longest. Write `.mp4` instead of `.gif` for smoother colours (needs ffmpeg).

The choice of $\kappa$ matters. With the default $\kappa = 0.005$, diffusion smooths the blob out within a few time units. A smaller $\kappa$ keeps the filaments sharp but needs a finer grid: keep the grid Péclet number printed by the solver below about 2.

## Larger runs on a bigger machine

Halving $\kappa$ with each doubling of the grid keeps the grid Péclet number at about 1.8 and makes the mixing sharper. The work then grows about 8× per doubling:

| Grid | `--kappa` | Steps per flow period ($t = 10$) | Cell updates per period | One core at 2.5e8 updates/s | Full-size snapshot |
| --- | --- | --- | --- | --- | --- |
| 512 × 512 | 0.001 | 21,150 | 5.5e9 | 0.4 min | 2 MB |
| 1024 × 1024 | 0.0005 | 42,300 | 4.4e10 | 3 min | 8 MB |
| 2048 × 2048 | 0.00025 | 84,600 | 3.6e11 | 24 min | 34 MB |
| 4096 × 4096 | 0.000125 | 169,200 | 2.8e12 | 3.2 h | 134 MB |

`heat_parallel` divides the time by roughly the number of cores times the parallel efficiency. On a laptop, 1024² ran 4.2× faster than serial with 2 ranks × 2 threads; measure the other machine with a short run (`--steps 1000`) before planning a long one.

[hpc/big_run.sh](hpc/big_run.sh) runs one such case with snapshots, serial or parallel; its settings come from environment variables:

```bash
# 2048^2 for two flow periods on 16 cores, snapshots as 512^2 block averages (8 MB each)
HEAT2D_REPO=$PWD RANKS=16 THREADS=1 nohup bash hpc/big_run.sh > outputs/run2048.log 2>&1 &
```

- `--snapshot-coarsen K` (`COARSEN`) stores K × K block averages, which keep the total mass; K must divide `nx` and `ny`.
- The final field goes to `final.npy`; binary output keeps large grids small (a 4096² CSV would be about 400 MB).
- `nohup` (or `tmux`) keeps the run going after logging out; the report is at the end of the log file.
- Animate on the big machine, or copy the snapshot folder back, for example `rsync -av machine:parallel_heat2d/outputs/run2048 outputs/`.
- Memory is not a constraint: about 5 × 8 × nx × ny bytes in total, 0.7 GB at 4096².

## Repository layout

| Path | Contents |
| --- | --- |
| [src/](src/README.md) | C++ sources: shared components (`common/`), serial solver (`serial/`), MPI + OpenMP solver (`parallel/`). |
| [tests/](tests/) | C++ unit tests, run by `ctest`. |
| [scripts/](scripts/README.md) | Correctness checks, plotting and animation; planned benchmark tools. |
| [config/defaults.json](config/defaults.json) | Default numerical and execution settings. |
| [experiments/](experiments/README.md) | Correctness and performance experiment configurations. |
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
- Measure strong scaling, weak scaling, rank/thread allocation and communication strategies.
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
