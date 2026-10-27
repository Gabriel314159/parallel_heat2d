# Verification record

Completed locally on 27 October 2026. These are measured software/numerical checks, not HPC performance conclusions.

## Environment

- Linux Ubuntu 24.04 container; GNU C++ 13.3.0.
- CMake 4.4.4; Release builds.
- MPICH 5.0.2; GCC OpenMP 4.5 support.
- Python 3.12; NumPy 2.3.5; Matplotlib 3.10.8.
- The restricted container needed MPICH's OFI/shared-memory transport and CMA disabled (`MPIR_CVAR_CH4_NETMOD=ofi`, `FI_PROVIDER=shm`, `MPIR_CVAR_CH4_CMA_ENABLE=0`). These testing-environment settings are not requirements for the documented Ubuntu/Open MPI setup.

## Completed checks

- Release compilation of the serial and MPI/OpenMP executables without compiler warnings.
- 10 checks in the main correctness tool, including seven serial/parallel full-field comparisons, a hand-computed first step, rejection of invalid numerical inputs and sine refinement.
- Maximum serial/parallel difference in the seven configured comparisons: 0 in this environment.
- 32 additional serial/parallel comparisons covering row and 2D decomposition, blocking and overlap, one/two threads, non-divisible dimensions, one-cell local dimensions and larger row/column messages.
- Maximum difference in the extended comparisons: 0.0.
- An impossible process grid was rejected without hanging.
- Independent sine-field refinement at final time 0.01, for grids 16, 32 and 64: maximum errors [0.00060343418872566, 0.00016828338462759884, 4.3884705987262684e-05]. The decreasing errors support numerical convergence; they are separate from serial/parallel equivalence tolerances.
- The small strong-scaling pipeline retained warmups/repetitions and generated runtime, speedup and efficiency PNGs.
- Smaller local weak-scaling and fixed-core strategy configurations exercised those runner and plotter branches. These short runs only check the workflow.
- A real temperature example was exported and plotted; see `examples/`.
- A separate serial-only Release build passed the serial numerical checks; dynamic-library inspection confirmed that it links neither MPI nor OpenMP.

## Scope still to verify externally

The Windows desktop/WSL installation and VS Code UI have not been executed in this container. Open MPI is the recommended Ubuntu installation path, while the actual MPI runtime tested here was MPICH. All communication uses standard MPI-3 APIs.

Atlas/Vanda queue configuration, multi-node networking, final resource binding and the full-size performance manifests still require runs in the actual HPC environment. The final report has not been completed. No claim of parallel speedup is made from tiny local checks.
