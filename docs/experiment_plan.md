# Validation and performance plan

Status: planned. No checks or benchmarks are claimed as completed.

## Correctness

### Serial and parallel equivalence

For each small case, compare the final full matrix from `heat_serial` with each selected parallel strategy. Use identical global grid dimensions, geometry, diffusivity, flow parameters, initial condition, time step, step count and precision. Compare with an absolute and relative tolerance, initially `1e-10` each, and record the observed error.

Test one rank, several ranks, multiple threads, non-divisible dimensions and very small local blocks. Include both row and 2D decompositions, and both blocking and overlap modes. The cases cover MPI only (several ranks, one thread), OpenMP only (one rank, several threads) and MPI + OpenMP, so that every layout of the performance experiments is checked before it is timed; `scripts/check_correctness.py` labels each run with its layout. Run these with the computed flow switched on (a short spin-up, then the pigment), mostly with the steady wind plus one case with the periodic wind.

### Mass conservation

With no-flux walls the total scalar $M = \Delta x\,\Delta y\sum c_{i,j}$ is conserved exactly by the flux-form scheme. For every correctness run, compute $M$ before and after the time loop and require

$$
\frac{\lvert M_{\text{end}} - M_0\rvert}{\Delta x\,\Delta y\sum_{i,j}\lvert c^0_{i,j}\rvert} \le 10^{-12}.
$$

For non-negative fields such as `pigment_drop` the denominator is $M_0$ itself, so this is the relative mass change; it stays meaningful for `cosine_decay`, where $M_0 = 0$. The solvers report this quantity as `mass_error`. A larger drift points to inconsistent face fluxes between neighbouring cells or ranks, a wrong ghost fill, or a velocity field that is not discretely divergence-free.

### Constant field

Initialize every cell to the same constant $c_0$ with the computed flow switched on. Because the discrete flow is divergence-free and the walls are no-flux, every update must retain $c_0$ up to floating-point roundoff. This checks initialization, the ghost fill, halo updates, the face velocities and the update formula together.

### Analytical cosine decay

With the flow switched off ($A = 0$), initialize

$$
c(x,y,0) = A_0\cos\!\Big(\frac{\pi x}{l_x}\Big)\cos\!\Big(\frac{\pi y}{l_y}\Big).
$$

Two separate comparisons:

1. **Exact discrete solution.** The cosine mode is an eigenvector of the discrete scheme with no-flux walls (see `docs/design.md`), so after $n$ steps the computed field must equal $(1-\lambda\Delta t)^n\,c(x,y,0)$ to rounding, with $\lambda$ given in the design. This is an implementation check and uses a tight relative tolerance, initially `1e-12`.
2. **Continuous solution.** Compare with $c(x,y,0)\,\exp\!\big(-\kappa\pi^2(1/l_x^2+1/l_y^2)\,t\big)$ at the actual $t = n\,\Delta t$. This difference is discretization error, not a bug, and is reported separately. Run a mesh-refinement study at a common physical final time, recomputing stable time steps, and check that the error decreases at the expected order. Record both time and spatial resolution.

The correctness manifest includes one cosine case for baseline equivalence and analytical error reporting. It does not by itself define the complete refinement study; add the common-final-time resolution series when implementing that study.

### Computed flow

These checks cover the wind-driven flow of `design.md`. Each one targets one part of the solver. They run in `ctest` (`tests/test_poisson.cpp`, `tests/test_flow.cpp`, `tests/test_serial_solver.cpp`) and in `scripts/check_correctness.py`.

**Rest stays at rest.** Without wind ($F_0 = 0$) and starting from $\zeta = 0$, the fields $\psi$ and $\zeta$ must stay exactly 0. This checks that no term creates flow out of nothing.

**Poisson solver.** With $\psi = 0$ on the walls,

$$
\psi(x,y) = \sin\!\Big(\frac{\pi x}{l_x}\Big)\sin\!\Big(\frac{\pi y}{l_y}\Big)
$$

is an exact eigenfunction of the discrete five-point Laplacian on the corners, so the CG solver must reproduce it to its tolerance. On a random right-hand side, the true relative residual of the solution must meet the tolerance. On a smooth right-hand side with a known solution, the error must fall by about 4× each time the grid is refined by 2×. A zero right-hand side must give $\psi = 0$ exactly, a warm start from the solution must need no iterations, and a solve that cannot reach the tolerance must stop with a message.

**Velocities from ψ.** Sampling an analytic streamfunction on the corners and applying the face formulas must give exactly the face velocities computed from the streamfunction directly. The discrete divergence must stay at rounding level for any $\psi$.

**Conservation of the Arakawa Jacobian.** Without wind and viscosity, the discrete energy $E$ and enstrophy $Z$ (`design.md`, "Diagnostics") must be conserved up to the time-stepping error, which must shrink at least like $\Delta t^2$ as $\Delta t$ is reduced (the forward-Euler start-up step is first order).

**Linear steady state.** Without the advection term and with the steady wind, the steady state solves $\nu\nabla^2\zeta = -F$, $\zeta = -\nabla^2\psi$, with $\psi = \zeta = 0$ on the walls. Expanding $1 = \sum_{m\ \text{odd}} \frac{4}{m\pi}\sin(m\pi x/l_x)$ gives the exact solution

$$
\psi(x,y) = \sum_{m\ \text{odd}} \frac{4F_0}{m\pi\,\nu\,\lambda_m^2}\,\sin\!\Big(\frac{m\pi x}{l_x}\Big)\sin\!\Big(\frac{2\pi y}{l_y}\Big),
\qquad
\lambda_m = \Big(\frac{m\pi}{l_x}\Big)^2 + \Big(\frac{2\pi}{l_y}\Big)^2 .
$$

The discrete steady state (two Poisson solves) must converge to it at second order; time-stepping from rest must reach the discrete steady state; and for $F_0 > 0$ the southern gyre must turn counterclockwise and the northern one clockwise.

**Manufactured solution.** We choose a smooth $\psi(x,y,t)$ that vanishes, with its Laplacian, on the walls, compute the source term that makes it an exact solution, and add that source to $F$. The error must converge at second order in space. This one check exercises every term at once.

**The two wind forms.** `--wind periodic` with $\alpha = 0$ must give exactly the same result as `--wind steady`.

**Pigment.** The pigment mass must stay conserved to about $10^{-13}$, and a constant pigment field must stay constant under the flow.

**Time-step check.** A run whose fixed $\Delta t$ is too large for the speeds the flow reaches must stop with a message naming the violated limit.

**Serial and parallel.** The CG solver's global sums add the rank contributions in a different order in parallel, so the serial and parallel solutions are no longer identical bit for bit. They must agree within `serial_parallel_atol` and `serial_parallel_rtol` (initially $10^{-10}$); on the correctness cases the largest difference is about $5\times10^{-14}$. Over long runs of an unsteady flow the two may drift apart, even though both are correct; long runs are then compared through statistics such as the time-mean kinetic energy.

## Timing protocol

- Compile all implementations in Release mode with comparable optimization flags and arithmetic.
- Run the serial reference and parallel variants on comparable compute hardware, ideally the same node type.
- Measure the time-stepping loop, including required MPI communication and synchronization.
- Keep allocation, initial setup, field gathering, output, plotting and checking outside the loop timer.
- Put an MPI barrier immediately before timing, then reduce rank-local elapsed times with `MPI_MAX` after the timed loop.
- Perform one warmup and five measured repetitions initially. Preserve every repetition and report a median plus variation.
- Record actual MPI ranks, OpenMP threads, node count, binding, libraries, hardware and Git commit. The solvers report the ranks, threads and nodes they used; the binding comes from the launcher and must be recorded separately.
- Allocate enough CPUs and memory for the largest layout. Avoid oversubscription and unrelated concurrent runs.

Optional phase timers should clearly distinguish wait time from complete communication time. Under overlap, diagnostic phase times are not additive estimates of total runtime.

## Performance experiments

The five experiments are listed in [experiments/README.md](../experiments/README.md). Each one changes one factor and keeps the others at the baseline (2D decomposition, blocking exchange and one thread per rank, unless the experiment changes it), so that a difference in time can be attributed to that factor. Experiments 1 to 4 solve the same fixed problem $P$ (`experiments/benchmark_problem.json`), so one serial reference serves all of them. They run only after the correctness checks have passed with the same build on the same machine.

The grid and the number of steps of $P$ are provisional. They are chosen so that a serial run takes a few minutes on a laptop (the first step took 1.7 s). With the conjugate gradient solver, the cost of a run grows about 8× each time the grid doubles in both directions, so time one run on Atlas or Vanda and adjust them in `benchmark_problem.json`.

```text
speedup S(p) = median_serial_time / median_parallel_time(p)
efficiency E(p) = S(p) / p
p = MPI ranks * OpenMP threads per rank
```

The denominator counts physical CPU cores when one thread is bound to each physical core. Record SMT use explicitly if tested. Do not describe a one-rank MPI run as the pure serial implementation.

### 1. Serial vs MPI (strong scaling over ranks)

Hold $P$ fixed and increase the number of ranks, with one thread per rank. Report the runtime, $S(p)$ and $E(p)$ against the pure serial program. The one-rank run shows what the MPI setup and the halo code cost on their own, when compared with the serial run.

### 2. Serial vs OpenMP (strong scaling over threads)

Hold $P$ fixed on one rank and increase the number of threads. Report the same quantities. The conjugate gradient loops read and write large arrays with little arithmetic per value, so memory bandwidth may limit the speedup. Bind each thread to a core and record the binding.

### 3. Pure MPI vs MPI + OpenMP

Fix the number of cores at 16 and change how they are split between ranks and threads, from 16 ranks × 1 thread to 1 rank × 16 threads. Report the runtime of each layout and its ratio to pure MPI. Inspect the local block dimensions and the exposed boundary lengths. A layout with less communication is not guaranteed to be faster if it changes message overhead, local threading efficiency or memory placement.

### 4. Communication strategy

At 16 ranks with one thread each, compare row and 2D decomposition, and blocking and overlapping exchange, in four runs. Rows give each rank two neighbours and long messages; 2D blocks give up to four neighbours and shorter messages, and the columns are packed before sending.

What this experiment can show is limited. `--exchange` changes only the halo exchange of the pigment, which happens once per step. The Poisson solve exchanges halos once per conjugate gradient iteration, hundreds of times per step, and always blocks. It takes about 99% of the run time (informal laptop measurement, 4 ranks, 256 × 512 cells: 99.2%), and the blocking and overlapping runs differed by less than the run-to-run noise. Expect no visible difference between blocking and overlap until the Poisson solve overlaps its own exchange. The decomposition, in contrast, changes the messages of every iteration. Report what is measured, and do not claim that overlap is faster before measuring it.

### 5. Weak scaling

Keep a 512-by-512 interior block and four threads per rank, then use 1, 4 and 16 ranks with square process grids. Keep `dx` and `dy` fixed and grow the physical domain according to the configuration rules. Keep the flow settings and the speed estimate `--max-speed` fixed, so that `dt` and the step count are the same for every run. The pigment and vorticity updates then do the same work per core in every run, but the Poisson solve does not: the number of CG iterations per step grows with the global grid (about doubling when the grid doubles in each direction). Record the CG iterations and the Poisson time, and report the time per CG iteration alongside the total time.

```text
weak-scaling efficiency = median_time(reference_resources) / median_time(larger_resources)
```

The reference here is one MPI rank with four threads, not the single-threaded serial executable. Report this baseline explicitly. Weak scaling does not compare the same physical domain at different resolutions.

## Planned figures

- Pigment field for a small demonstration case, with an optional velocity overlay.
- Maximum absolute serial/parallel difference for validation cases.
- Relative mass change for the validation cases.
- Analytical error versus resolution for the cosine case.
- Runtime, speedup and efficiency versus ranks (experiment 1) and versus threads (experiment 2).
- Runtime of each rank/thread layout at 16 cores (experiment 3) and of the four decomposition and exchange variants (experiment 4).
- Weak-scaling efficiency and time per CG iteration versus ranks (experiment 5).

Use only retained, measured data. Label the machine, grid, steps, resource definition and timing scope. Separate observations from proposed explanations.
