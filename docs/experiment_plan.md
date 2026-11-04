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
- Bind one thread to each physical core, and leave the other hardware threads of every core unused (no SMT or hyper-threading; `hpc/README.md`, "Binding"). A run counts only if the solver reports its `binding` as `cores`.
- Record actual MPI ranks, OpenMP threads, node count, binding, libraries, hardware and Git commit. The solvers report the ranks, threads and nodes they used, the fewest physical cores any rank may run on (`cores_per_rank`) and the `binding`; `heat_parallel` warns when a rank has fewer cores than threads or when the binding is not `cores`. Record the launcher's binding options as well.
- Allocate enough cores and memory for the largest layout. Avoid oversubscription and unrelated concurrent runs.
- On Vanda a node has 72 cores in two sockets of 36 (`hpc/README.md`). Every timed run reserves its nodes exclusively (`place=excl`); a run on several nodes uses 64 cores of each, 32 per socket, so that 256 cores are four nodes and every layout of experiment 3 fits the sockets. Mark in the figures where the runs start to span nodes.

The parallel solver reports two phase times of the Poisson solves: `poisson_halo_s`, the time inside the CG's halo exchanges (under overlap only posting and waiting, the part not hidden), and `poisson_sum_s`, the time inside its global sums, including the wait for the slowest rank. Each is the largest over the ranks, so they are diagnostic and not additive estimates of the total runtime.

## The whole computation and its dominant part

We discuss every scaling result at two levels: first the whole computation, then the parts that dominate it, as the feedback on the progress report asked. The solvers report what we need:

- $T$, the time of the whole time-stepping loop (`elapsed_s`);
- $T_P$, the part of it spent in the Poisson solves (`poisson_s`);
- $n$, the number of steps, and $\bar k$, the mean number of CG iterations per step (`poisson_iterations_mean`);
- inside $T_P$, the time in the CG's halo exchanges and in its global sums (`poisson_halo_s`, `poisson_sum_s`).

1. **The whole computation.** The runtime, speedup and efficiency of every experiment are computed from $T$. These are the main results.
2. **The dominant part.** The Poisson solve takes about 99% of $T$ on the laptop, so we analyse it separately. At every core count we report its share $T_P/T$, the iterations per step $\bar k$, and the time per CG iteration

   $$
   t_{\text{it}} = \frac{T_P}{n\,\bar k}.
   $$

   $t_{\text{it}}$ separates the parallel implementation (the cost of one iteration) from the algorithm (more iterations on a bigger grid). The halo and sum times then show what limits $t_{\text{it}}$ at many cores: the computation, the messages between neighbours, or the global sums.
3. **The rest.** The pigment update, the vorticity update and the velocities take $T - T_P$. We report its share as well. If it grows to a sizeable part of $T$ at many cores, we analyse it in the same way.

All these times are maxima over the ranks, so the parts do not add up exactly, and the shares are estimates.

## Performance experiments

The five experiments are listed in [experiments/README.md](../experiments/README.md). Each one changes one factor and keeps the others at the baseline (2D decomposition, blocking exchange and one thread per rank, unless the experiment changes it), so that a difference in time can be attributed to that factor. Experiments 1 to 4 solve the same fixed problem $P$ (`experiments/benchmark_problem.json`), so one serial reference serves all of them. They run only after the correctness checks have passed with the same build on the same machine.

$P$ is the $1 \times 1$ basin on $1024 \times 1024$ cells; its number of steps is provisional. A serial run of 100 steps takes about 8 minutes on a laptop (about 840 CG iterations per step from rest, at about 5.5 ns per cell and iteration). The flow starts at rest and hardly moves in 100 steps, so the iterations per step are those of a starting flow; the report says so. A developed flow would need a restart from a saved state, which the solvers do not have. With the conjugate gradient solver, the cost of a run grows about 8× each time the grid doubles in both directions, so time one run on Vanda and adjust the steps in `benchmark_problem.json`.

```text
speedup S(p) = median_serial_time / median_parallel_time(p)
efficiency E(p) = S(p) / p
p = MPI ranks * OpenMP threads per rank
```

The number of cores $p$ counts physical cores: every thread is bound to a core of its own, and no core runs two of our threads. Do not describe a one-rank MPI run as the pure serial implementation.

**Two serial references.** On the laptop (Apple clang), `heat_serial` took twice as long per CG iteration as `heat_parallel` on one rank and one thread: 0.98 s against 0.49 s of Poisson time on $256 \times 512$ cells, with the same 309 iterations per step. The same serial solver compiled with the OpenMP pragmas and run on one thread took 0.49 s, so the gap comes from how the compiler optimises the shared CG loops without OpenMP, not from MPI. We therefore keep both: `heat_serial`, the pure serial program, and `heat_serial_omp`, the same solver compiled like `heat_parallel` and forced to one thread. Every speedup is reported against both, $S(p)$ against `heat_serial` and $S_{\text{omp}}(p)$ against `heat_serial_omp`; the second separates the parallel gain from the code generation. Check on Vanda with GCC whether the gap is there at all, and report it either way.

The course asks for scaling tests up to at least 256 cores. Experiments 1 and 5 go that far, and experiments 3 and 4 compare their variants on 256 cores; experiment 2 stops at one node, the limit of OpenMP. On Vanda, 256 cores are four nodes, in the charged queue `batch_cpu` (the free queue stops at 160 cores).

### 1. Serial vs MPI (strong scaling over ranks)

Hold $P$ fixed and increase the number of ranks from 1 to 256 in powers of 2, with one thread per rank. Report the runtime, $S(p)$ and $E(p)$ against the pure serial program. The one-rank run shows what the MPI setup and the halo code cost on their own, when compared with the serial run. Report `poisson_halo_s` and `poisson_sum_s` next to the runtime: they show how much of the Poisson time goes to communication as ranks are added, and so where the speedup stops growing.

At 256 ranks a rank holds only $64 \times 64$ cells of $P$, but every CG iteration still needs two global sums over all ranks, and their cost grows with the number of ranks. We expect the speedup to stop growing where the time of the sums reaches the time of the computation. If that happens early, a larger $P$ would show more of the curve, at the price of a longer serial reference (about 8× per doubling of the grid in both directions; fewer steps would compensate). Experiments 3 and 4 also run $P$ on 256 cores and would change with it. Decide after timing $P$ on the cluster.

### 2. Serial vs OpenMP (strong scaling over threads)

Hold $P$ fixed on one rank and increase the number of threads from 1 to the 72 cores of a Vanda node: 1, 2, 4, 8, 16, 32, 36, 64 and 72. OpenMP threads share the memory of one node, so this experiment cannot go beyond one node; experiments 1 and 5 go to 256 cores.

Report the same quantities. The conjugate gradient loops read and write large arrays with little arithmetic per value, so memory bandwidth may limit the speedup. Bind each thread to a core of its own (`binding` `cores`): a launcher that binds the single rank to one core would otherwise run all threads there. With close binding the first 36 threads fill one socket; mark 36, where the threads start to use the second one. Memory placement matters here, because Linux puts a page on the socket of the thread that writes it first. The first sweep (`strong_scaling_openmp.csv`) zeroed every field from the main thread: the time stopped falling at 36 threads (10.7 s) and was 11.4 s at 64 and 11.8 s at 72, which we expect comes from the threads of the second socket reading memory next to the first. `Field` now zeroes its storage in parallel, with the schedule of the solver loops (`docs/design.md`, "Data storage and shared computation"), and the sweep was repeated (`strong_scaling_openmp_firsttouch.csv`). We report both, as the cost of the memory placement.

### 3. Pure MPI vs MPI + OpenMP

Fix the number of cores at 256 and change how they are split between ranks and threads, from 256 ranks × 1 thread (pure MPI) to 4 ranks × 64 threads, in powers of 2. The threads of a rank share its memory, so they must fit on one node. Every layout runs on four Vanda nodes with 64 cores each, 32 per socket: 8 × 32 is one rank per socket and 4 × 64 one rank per node. Pure OpenMP is experiment 2.

| Ranks × threads | Process grid | Cells per rank |
| --- | --- | --- |
| 256 × 1 | 16 × 16 | 64 × 64 |
| 128 × 2 | 16 × 8 | 128 × 64 |
| 64 × 4 | 8 × 8 | 128 × 128 |
| 32 × 8 | 8 × 4 | 256 × 128 |
| 16 × 16 | 4 × 4 | 256 × 256 |
| 8 × 32 | 4 × 2 | 512 × 256 |
| 4 × 64 | 2 × 2 | 512 × 512 |

Report the runtime of each layout and its ratio to pure MPI, then the Poisson quantities of "The whole computation and its dominant part". The edges of a rank's block set the length of its messages. The MPI calls run between the OpenMP loops, on one thread of each rank, while the others wait. With fewer ranks there are fewer and longer messages, and fewer ranks in each global sum. In exchange, one thread communicates for the whole rank, each of the three OpenMP loops of a CG iteration ends by waiting for all threads of the rank, and a rank whose threads span two sockets (4 × 64) depends on where its memory sits (experiment 2). That layout was run again after the first-touch change (`hybrid_layouts_firsttouch.csv`); the first run stays as the "before". A layout with less communication is therefore not guaranteed to be faster.

**On GPUs.** We also run $P$ with its Poisson solve on GPUs, through OpenMP target offload, with one rank per GPU: 1 rank × 36 threads with one GPU, then 2 ranks × 36 threads with the two GPUs of one node. A Vanda GPU node has two NVIDIA A40, and the site gives each GPU one socket of 36 cores, on which the rest of the step runs. Compare the runtime with the CPU layouts on 256 cores and with one CPU node (72 cores), and report the time of the copies between host and GPU. The A40 computes in double precision at 1/64 of its single-precision rate, but the CG moves much more data than it computes, so its memory bandwidth (about 700 GB/s) should limit it first. The two-GPU run needs the charged queue `batch_gpu` (the free queue allows 36 cores). On an A40, `nvc++` (NVHPC 26.1) and `clang++` (LLVM 20.1) both ran an OpenMP target region on the GPU; `g++` cannot offload there (`hpc/README.md`). The GPU build is described in `docs/design.md` ("The Poisson solve on a GPU"); it uses the blocking exchange only, and reports `poisson_device` and the copy time `poisson_copy_s`. Check `poisson_device` = `gpu` in every GPU run. Run the blocking correctness cases with the GPU build on a GPU node before timing it.

### 4. Communication strategy

At 256 ranks with one thread each, compare row and 2D decomposition, and blocking and overlapping exchange, in four runs. The runs span several nodes, so part of the messages cross the network. Rows give each rank two neighbours and long messages; 2D blocks give up to four neighbours and shorter messages, and the columns are packed before sending. On 256 ranks the messages differ a lot. With rows, a rank of $P$ holds 4 rows of 1024 cells and exchanges 1024 values with each of its two neighbours; a 2D block of $64 \times 64$ cells exchanges 64 values with each of its four. Rows thus send eight times as many values, in half as many messages.

`--exchange` applies to both halo exchanges: the pigment's, once per step, and the conjugate gradient's, once per iteration, hundreds of times per step (`docs/design.md`, "Halo exchange"). The Poisson solve takes about 99% of the run time, so the CG exchange is the one that matters. The first version overlapped only the pigment exchange, and blocking and overlapping runs then differed by less than the run-to-run noise (informal laptop measurement, 4 ranks, 256 × 512 cells).

Both modes post the same messages at once; blocking waits for them before computing, overlap computes the inner corners first. The two global sums of every CG iteration always wait, so overlap can save at most the CG halo time of the blocking runs (`poisson_halo_s`). Compare that with the difference in runtime. While the messages travel, overlap computes the corners that read no process halo: about 94% of a $64 \times 64$ block, but only half of a block of 4 rows. Within one node the messages are copies through shared memory, and the gain may be small: an informal laptop test with the CG overlap (4 ranks, 256 × 512 cells) showed halo exchanges taking 6–11% of the Poisson time and no difference beyond the noise. On 256 ranks part of the messages cross the network, where there may be more to hide. The decomposition changes the messages of every iteration. Report what is measured, and do not claim that overlap is faster before measuring it.

### 5. Weak scaling

We start from one core and give every added core the same work: one rank per core, one thread per rank, and a block of $256 \times 256$ cells per rank. The process grid is square, $q \times q$ with $q = 1, 2, 4, 8, 16$, so the runs use

$$
p = q^2 = 1,\ 4,\ 16,\ 64,\ 256 \ \text{cores}
$$

on a global grid of $256q \times 256q$ cells. We keep $\Delta x = \Delta y = 1/256$ fixed, so the basin grows to $q \times q$ and keeps the $1 \times 1$ shape of the benchmark. The run on one core has the default grid ($256 \times 256$ cells on the $1 \times 1$ basin), with the settings of $P$ otherwise; both serial references run it as well, which shows what MPI costs on one core. The flow settings and the speed estimate `--max-speed` stay fixed, so $\Delta t$ and the number of steps are the same in every run.

We keep one thread per rank, so that the number of cores is the number of ranks and the series starts on a single core. How to split the cores between ranks and threads is the question of experiment 3.

The pigment and vorticity updates then do the same work per core in every run, but the Poisson solve does not: the CG iterations per step grow with the global grid, about doubling each time $q$ doubles. An informal laptop test (10 steps from rest) gave 293 iterations per step at $q = 1$ and 542 at $q = 2$, so we expect about 4,000 at $q = 16$. The total time therefore grows even if the parallel implementation scales perfectly. We report two efficiencies,

$$
E_w(p) = \frac{T(1)}{T(p)},
\qquad
E_{w,\text{it}}(p) = \frac{t_{\text{it}}(1)}{t_{\text{it}}(p)},
$$

where $T$ is the median loop time and $t_{\text{it}}$ the time per CG iteration (see "The whole computation and its dominant part"). $E_w$ is the weak scaling of the whole computation. $E_{w,\text{it}}$ removes the growth of the iterations, so it measures how well one CG iteration scales. The reference $T(1)$ is the one-rank run of `heat_parallel`; report it next to the serial run.

At 256 cores a run of 100 steps then needs about 400,000 CG iterations. At about 0.3 ms per iteration (the laptop's 0.6 ms for 4 ranks with twice the cells per rank), that is about 2 minutes per run. Time the 64-core run first; if the runs are too long, reduce the steps for this experiment (a `problem` object with `steps` in `weak_scaling.json`). Weak scaling does not compare the same physical domain at different resolutions.

## Planned figures

- Pigment field for a small demonstration case, with an optional velocity overlay.
- Maximum absolute serial/parallel difference for validation cases.
- Relative mass change for the validation cases.
- Analytical error versus resolution for the cosine case.
- Runtime, speedup and efficiency versus ranks (experiment 1, up to 256) and versus threads (experiment 2, up to one node), with the node and socket boundaries marked.
- Time inside the CG's halo exchanges and global sums versus ranks (experiment 1).
- Share of the Poisson solve, $T_P/T$, and time per CG iteration, $t_{\text{it}}$, versus cores (experiments 1, 2 and 5).
- Runtime of each rank/thread layout on 256 cores and of the runs with one and two GPUs (experiment 3), and runtime and CG halo time of the four decomposition and exchange variants (experiment 4).
- Both weak-scaling efficiencies, $E_w$ and $E_{w,\text{it}}$, and the CG iterations per step versus cores (experiment 5).

Use only retained, measured data. Label the machine, grid, steps, resource definition and timing scope. Separate observations from proposed explanations.
