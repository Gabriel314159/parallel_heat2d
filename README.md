# Parallel Heat2D

A course project on a parallel two-dimensional advection–diffusion solver (pigment in a double-gyre flow) using C++, MPI and OpenMP.

We will implement the same numerical method in serial and in parallel, then compare how domain decomposition, halo communication and thread allocation affect correctness and performance.

## Current status

**This is a repository scaffold. The numerical solvers and experiment tools have not been implemented yet.**

Included now:

- Repository description, project scope and development conventions.
- CMake configuration with solver builds disabled by default.
- Default parameters and planned experiment configurations.
- Design notes, an experiment plan and a report outline.
- Source-directory guides and an HPC job template.

Not included yet:

- C++ solver source files or executable binaries.
- Benchmark runners, correctness-checking scripts or plotting scripts.
- Measured performance results or completed validation results.

The solver commands below describe the **planned interface**. They cannot run until the relevant source files are implemented.

For initial setup instructions in Chinese, see [START_HERE_zh.md](START_HERE_zh.md).

Initial repository/configuration checks are recorded in [docs/scaffold_checks.md](docs/scaffold_checks.md). These are not numerical solver tests.

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

## Planned implementations

| Executable or mode | Purpose |
| --- | --- |
| `heat_serial` | A single-threaded reference implementation with no MPI dependency. |
| `heat_parallel --decomp rows --exchange blocking` | MPI row decomposition with optional OpenMP computation. |
| `heat_parallel --decomp 2d --exchange blocking` | A two-dimensional MPI process grid with optional OpenMP computation. |
| `heat_parallel --decomp rows/2d --exchange overlap` | Nonblocking halo exchange with interior computation before waiting. |

The actual `--decomp` values will be `rows` or `2d`; `rows/2d` above denotes two choices. Set `OMP_NUM_THREADS=1` for an MPI-only baseline. A single-rank run of `heat_parallel` still includes MPI setup and is not the pure serial reference.

## Repository layout

| Path | Contents |
| --- | --- |
| [src/](src/README.md) | Planned C++ sources: shared components (`common/`), serial solver (`serial/`) and MPI + OpenMP solver (`parallel/`). |
| [config/defaults.json](config/defaults.json) | Default numerical and execution settings. |
| [experiments/](experiments/README.md) | Planned correctness and performance experiment configurations. |
| [scripts/](scripts/README.md) | Specifications for the tools to implement later. |
| [hpc/](hpc/README.md) | HPC setup notes and an unconfigured job template. |
| [results/raw/](results/raw/README.md) | Measured benchmark CSV files, when available. |
| [results/figures/](results/figures/README.md) | Selected figures generated from measured data. |
| [docs/design.md](docs/design.md) | Numerical and software design decisions. |
| [docs/experiment_plan.md](docs/experiment_plan.md) | Correctness checks and fair timing conventions. |
| [docs/ai_usage.md](docs/ai_usage.md) | Record of AI assistance and subsequent verification. |
| [docs/references.md](docs/references.md) | Technical references. |
| [report/report.md](report/report.md) | Report outline; no experimental findings are claimed. |

## Requirements

- CMake 3.16 or newer.
- A C++17 compiler.
- An MPI implementation for the parallel target.
- Compiler support for OpenMP when OpenMP is enabled.
- Python 3.9 or newer for future experiment utilities.

The C++ solver will not depend on Python. `requirements.txt` lists NumPy and Matplotlib for future checking and plotting tools; these tools are not included yet. Configurations are JSON for the planned Python runner. The C++ programs will receive command-line arguments, so they do not need a JSON library.

## Configure the scaffold now

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

This checks the CMake project setup and compiler. It deliberately creates no solver executables because `HEAT2D_BUILD_SOLVERS` defaults to `OFF`. It does not check MPI or OpenMP availability in scaffold mode.

## Build after implementing the solver sources

**Future workflow:**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DHEAT2D_BUILD_SOLVERS=ON
cmake --build build --parallel
```

Enabling solver builds now will stop with a message listing the source files that still need to be implemented. MPI and OpenMP checks happen only after these files are present.

For a future serial-only build, also pass `-DHEAT2D_ENABLE_MPI=OFF`. A compiler wrapper can be supplied with `-DCMAKE_CXX_COMPILER=mpicxx` when appropriate for the environment.

## Planned command-line interface

Serial:

```bash
./build/heat_serial --nx 256 --ny 256 --steps 1000 --case pigment_drop --cfl 0.8 --output outputs/serial.csv
```

Parallel, with four MPI ranks and four threads per rank:

```bash
OMP_NUM_THREADS=4 mpiexec -n 4 ./build/heat_parallel --nx 256 --ny 256 --steps 1000 --case pigment_drop --cfl 0.8 --decomp 2d --exchange overlap --output outputs/parallel.csv
```

These examples use a POSIX shell. Native Windows shell syntax differs. HPC runs must use the scheduler and launcher approved for Atlas or Vanda, with matching CPU allocations.

Additional planned options are documented in [config/README.md](config/README.md). No benchmark runner command is provided yet because that tool does not exist.

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

1. Implement shared parameters, grids, initialization and the update formula.
2. Implement and validate `heat_serial`.
3. Implement row decomposition and blocking MPI halo exchange.
4. Add OpenMP computation and verify hybrid configurations.
5. Add two-dimensional decomposition, including column communication.
6. Add nonblocking exchange and computation overlap.
7. Implement the experiment utilities, run local checks and collect HPC measurements.
8. Write conclusions from the measured data and finish the report.

## Submission scope

The assignment asks for a PDF report, source code with running instructions, measured runtime figures from Atlas or Vanda, correctness validation and an explanation of AI use. The supplied assignment specifies submission before 13 November 2026 and mentions an optional progress report before 30 October. Confirm the precise cutoff time on Canvas.
