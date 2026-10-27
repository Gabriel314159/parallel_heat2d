# Parallel Heat2D

Two-dimensional heat diffusion in C++17, with a serial reference and MPI + OpenMP implementations. The same five-point finite-difference update is used by every solver.

## What is implemented

- Contiguous double-precision grids, fixed physical boundaries and double buffering.
- `heated_top`, `constant` and analytical `sine_decay` initial-condition cases.
- Automatic stable time-step selection and validation of explicit time steps.
- A pure serial executable, independent of MPI and OpenMP.
- MPI row or 2D block decomposition, including remainders and thin local blocks.
- Blocking halo exchange, and nonblocking exchange overlapping interior computation.
- OpenMP point updates with MPI calls confined to the initializing thread.
- Final-field CSV export, one-line JSON runtime/metadata output and full-field validation.
- Configuration-driven benchmark runs with warmups, raw repetitions, logs and provenance.
- Temperature, runtime, speedup, strong-efficiency and weak-efficiency plotting.
- VS Code build tasks and a serial GDB debug configuration for WSL Ubuntu.

**Local numerical checks have run. Atlas/Vanda measurements and the final report are still pending.** See [verification.md](docs/verification.md) for the actual tested environment and results. No local demo timing is presented as an HPC benchmark.

For Windows setup and detailed Chinese instructions, start with [START_HERE_zh.md](START_HERE_zh.md).

## Model

The plate has size `lx` by `ly`, constant positive diffusivity `alpha`, no source term, and fixed Dirichlet boundaries. The model is:

```text
dT/dt = alpha * (d2T/dx2 + d2T/dy2)

next(i,j) = current(i,j)
          + rx * (current(i,j-1) - 2*current(i,j) + current(i,j+1))
          + ry * (current(i-1,j) - 2*current(i,j) + current(i+1,j))
```

`i` is a top-down row index; `j` is a column index. `nx` and `ny` count interior columns and rows. Physical boundaries add two rows and two columns to a full field. Local MPI fields additionally hold one-cell halos; process halos are excluded from global output.

```text
dx = lx/(nx+1), dy = ly/(ny+1)
rx = alpha*dt/dx^2, ry = alpha*dt/dy^2
rx + ry <= 0.5

dt_max = 1 / (2*alpha*(1/dx^2 + 1/dy^2))
automatic dt = cfl * dt_max
```

Default `cfl=0.8` is a fraction of the maximum stable time step. It is not `rx` or `ry`. Both solvers derive spacing from global dimensions. Use two arrays: read current, write next, complete the step, then swap buffers.

The default demonstration uses 20°C initially, a 100°C top edge, and three 20°C edges. These are synthetic conditions, not measured material properties. Left/right edge values take precedence at the four corners.

## Build on Ubuntu / WSL

Install the system packages once:

```bash
sudo apt update
sudo apt install -y build-essential cmake gdb openmpi-bin libopenmpi-dev python3-venv python3-pip
```

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DHEAT2D_BUILD_SOLVERS=ON
cmake --build build --parallel 2
```

The explicit `HEAT2D_BUILD_SOLVERS=ON` overrides the old scaffold's cached `OFF` setting. MPI and OpenMP are enabled by default. To build only the serial executable, use a separate build directory with `-DHEAT2D_ENABLE_MPI=OFF`. To build MPI without OpenMP, pass `-DHEAT2D_ENABLE_OPENMP=OFF` and request one thread per rank.

CMake 3.16+, a C++17 compiler and an MPI-3+ implementation are sufficient for the core. Python 3.9+, NumPy and Matplotlib serve the experiment tools. The C++ solvers do not depend on Python or a JSON library.

## Run the solvers

Serial:

```bash
./build/heat_serial --nx 128 --ny 128 --steps 2000 --case heated_top --output outputs/serial.csv
```

Parallel on local Open MPI, two ranks and two threads per rank:

```bash
OMP_NUM_THREADS=2 mpiexec --bind-to none -n 2 ./build/heat_parallel --nx 128 --ny 128 --steps 2000 --case heated_top --decomp 2d --exchange overlap --output outputs/parallel.csv
```

The second example requests four CPU workers. Use `OMP_NUM_THREADS=1` for a two-worker demonstration. `--bind-to none` is an Open MPI option for this local correctness/demo command; configure deliberate placement for performance experiments.

Run `./build/heat_serial --help` for all options. `heat_parallel` accepts `--decomp rows|2d`, `--exchange blocking|overlap`, optional `--py`/`--px`, and an optional `--threads N` override. MPI ranks are set by the launcher. Without an explicit thread count or `OMP_NUM_THREADS`, the parallel executable uses one thread per rank.

Both commands print one JSON record. `elapsed_s` covers the solver loop only, including required communication and threading. The parallel record reports the maximum rank-local elapsed time. Initialization, halo datatype creation, gathering, checking and file output are outside that timer. Zero steps are supported for initialization checks.

## Check correctness

Fast local checks, at most two ranks and one thread per rank:

```bash
python scripts/check_correctness.py --quick
```

Full default checks, up to six ranks and two threads per rank:

```bash
python scripts/check_correctness.py
```

Extended MPI cases, including thin grids, all four strategies and larger row/column messages:

```bash
python tests/test_integration.py
```

For a serial-only build, use `python scripts/check_correctness.py --serial-only --build-dir build-serial`. Checks include a hand-computed first step, fixed boundaries, the discrete maximum principle, constant fields, full serial/parallel equivalence, invalid inputs and sine-field refinement at a common final time. Summaries and diagnostic fields are written under ignored `outputs/` directories.

Serial/parallel comparison tolerances are separate from analytical discretization errors. Parent directories for CSV outputs are created automatically. CSV matrices have numeric values only, with the full `(ny+2, nx+2)` physical field and round-trip double precision.

## View a temperature field

```bash
python scripts/plot_temperature.py outputs/serial.csv --output outputs/temperature.png --vmin 20 --vmax 100
```

Open the PNG in VS Code. [examples/](examples/README.md) contains a real local demonstration and a compact validation summary.

## Run experiments

Start with the small local pipeline check:

```bash
python scripts/run_benchmarks.py --config experiments/smoke_scaling.json
python scripts/plot_results.py outputs/smoke/metrics.csv --out-dir outputs/smoke/figures
```

Resolve a larger plan without running it:

```bash
python scripts/run_benchmarks.py --config experiments/strong_scaling.json --dry-run
```

The three full performance manifests cover strong scaling, weak scaling and fixed-core strategy comparisons. Run them only inside suitable CPU/memory allocations. `--allocated-cores N` declares the total available CPUs for a multi-node job; the default guard uses the local CPU affinity. This is a resource guard, not a scheduler allocation request.

Open MPI hybrid placement can be supplied per configuration:

```bash
OMP_PLACES=cores OMP_PROC_BIND=close python scripts/run_benchmarks.py --config experiments/strategies.json --mpi-arg=--map-by '--mpi-arg=slot:PE={threads}' --mpi-arg=--bind-to --mpi-arg=core
```

`{threads}` and `{ranks}` in launcher arguments are expanded for each run. Placement syntax is MPI/platform-specific. The same parser accepts `--launcher PATH` and repeated `--mpi-arg=VALUE`; a scheduler-based launcher must be configured for the actual HPC environment.

Raw CSV records distinguish warmups from measured repetitions. The plotter defaults to the latest session, reports measured medians and, for runtime plots, min–max ranges. Strong-scaling speedup uses the matching pure serial reference. Weak-scaling efficiency uses the smallest MPI layout with the same thread count, not the serial executable.

See [experiment_plan.md](docs/experiment_plan.md), [experiment configurations](experiments/README.md) and [HPC setup](hpc/README.md).

## Code map

| Path | Responsibility |
| --- | --- |
| `src/common/` | CLI, grid, initial/boundary conditions, stencil, CSV and JSON output. |
| `src/serial/` | Serial main program and time loop. |
| `src/parallel/` | MPI topology, halo exchange, threaded kernels, gather and parallel driver. |
| `scripts/common.py` | Configuration merging, launching, result parsing and metadata. |
| `scripts/check_correctness.py` | Numerical validation and refinement. |
| `scripts/run_benchmarks.py` | Sequential experiment runs and retained records. |
| `scripts/plot_temperature.py`, `scripts/plot_results.py` | Headless PNG plotting. |
| `.vscode/` | WSL build, validation and serial-debug tasks. |
| `tests/` | Extended MPI integration checks. |
| `hpc/` | Platform-adaptation notes and a job template. |
| `docs/`, `report/` | Design, verification, AI record and report outline. |

The 2D process-grid selector minimizes estimated local perimeter among fitting factorizations. Cartesian rank reordering is disabled; all topology operations still use the Cartesian communicator. Columns use `MPI_Type_vector` in both exchange modes. MPI calls stay on the initializing thread with `MPI_THREAD_FUNNELED`. The overlap path completes all requests before swapping buffers and handles blocks without a strict interior.

## Current limits

This version supports CPU computation, constant diffusivity, fixed boundaries and fixed step counts. It writes final fields; intermediate snapshots remain disabled (`snapshot_interval=0`). Global CSV gathering uses MPI `int` counts and is limited to at most `INT_MAX` interior points, with enough root memory. Large performance runs should avoid field export.

The included HPC template still needs actual queue/module/launcher settings. Numerical checks have run in a local Linux container; the complete Windows/WSL desktop workflow and Atlas/Vanda jobs have not been executed here. The report outline still needs real HPC measurements and discussion. AI assistance is recorded in [docs/ai_usage.md](docs/ai_usage.md).
