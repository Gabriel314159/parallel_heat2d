# Validation and performance plan

Status: planned. No checks or benchmarks are claimed as completed.

## Correctness

### Serial and parallel equivalence

For each small case, compare the final full matrix from `heat_serial` with each selected parallel strategy. Use identical global grid dimensions, geometry, diffusivity, flow parameters, initial condition, time step, step count and precision. Compare with an absolute and relative tolerance, initially `1e-10` each, and record the observed error.

Test one rank, several ranks, multiple threads, non-divisible dimensions and very small local blocks. Include both row and 2D decompositions, and both blocking and overlap modes. Run these with the time-dependent double gyre switched on, plus one steady ($\varepsilon = 0$) case to exercise the fill-once velocity path.

### Mass conservation

With no-flux walls the total scalar $M = \Delta x\,\Delta y\sum c_{i,j}$ is conserved exactly by the flux-form scheme. For every correctness run, compute $M$ before and after the time loop and require

$$
\frac{\lvert M_{\text{end}} - M_0\rvert}{\Delta x\,\Delta y\sum_{i,j}\lvert c^0_{i,j}\rvert} \le 10^{-12}.
$$

For non-negative fields such as `pigment_drop` the denominator is $M_0$ itself, so this is the relative mass change; it stays meaningful for `cosine_decay`, where $M_0 = 0$. The solvers report this quantity as `mass_error`. A larger drift points to inconsistent face fluxes between neighbouring cells or ranks, a wrong ghost fill, or a velocity field that is not discretely divergence-free.

### Constant field

Initialize every cell to the same constant $c_0$ with the double gyre switched on. Because the discrete flow is divergence-free and the walls are no-flux, every update must retain $c_0$ up to floating-point roundoff. This checks initialization, the ghost fill, halo updates, the face velocities and the update formula together.

### Analytical cosine decay

With the flow switched off ($A = 0$), initialize

$$
c(x,y,0) = A_0\cos\!\Big(\frac{\pi x}{l_x}\Big)\cos\!\Big(\frac{\pi y}{l_y}\Big).
$$

Two separate comparisons:

1. **Exact discrete solution.** The cosine mode is an eigenvector of the discrete scheme with no-flux walls (see `docs/design.md`), so after $n$ steps the computed field must equal $(1-\lambda\Delta t)^n\,c(x,y,0)$ to rounding, with $\lambda$ given in the design. This is an implementation check and uses a tight relative tolerance, initially `1e-12`.
2. **Continuous solution.** Compare with $c(x,y,0)\,\exp\!\big(-\kappa\pi^2(1/l_x^2+1/l_y^2)\,t\big)$ at the actual $t = n\,\Delta t$. This difference is discretization error, not a bug, and is reported separately. Run a mesh-refinement study at a common physical final time, recomputing stable time steps, and check that the error decreases at the expected order. Record both time and spatial resolution.

The correctness manifest includes one cosine case for baseline equivalence and analytical error reporting. It does not by itself define the complete refinement study; add the common-final-time resolution series when implementing that study.

## Timing protocol

- Compile all implementations in Release mode with comparable optimization flags and arithmetic.
- Run the serial reference and parallel variants on comparable compute hardware, ideally the same node type.
- Measure the time-stepping loop, including required MPI communication and synchronization.
- Keep allocation, initial setup, field gathering, output, plotting and checking outside the loop timer.
- Put an MPI barrier immediately before timing, then reduce rank-local elapsed times with `MPI_MAX` after the timed loop.
- Perform one warmup and five measured repetitions initially. Preserve every repetition and report a median plus variation.
- Record actual MPI ranks, OpenMP threads, node count, binding, libraries, hardware and Git commit.
- Allocate enough CPUs and memory for the largest layout. Avoid oversubscription and unrelated concurrent runs.

Optional phase timers should clearly distinguish wait time from complete communication time. Under overlap, diagnostic phase times are not additive estimates of total runtime.

## Strong scaling

Hold the full discrete problem fixed and increase resources. The initial manifest isolates rank scaling at one thread per rank. Measure a pure serial reference separately.

```text
speedup S(p) = median_serial_time / median_parallel_time(p)
efficiency E(p) = S(p) / p
p = MPI ranks * OpenMP threads per rank
```

The denominator counts physical CPU cores when one thread is bound to each physical core. Record SMT use explicitly if tested. Do not describe a one-rank MPI run as the pure serial implementation.

## Weak scaling

Keep a 512-by-512 interior block and four threads per rank, then use 1, 4 and 16 ranks with square process grids. Keep `dx` and `dy` fixed and grow the physical domain according to the configuration rules. Scale the double-gyre amplitude with the domain size so that the speed bounds $u_{\max}$, $v_{\max}$ stay fixed. Thus the automatic `dt`, step count and work per core remain comparable.

```text
weak-scaling efficiency = median_time(reference_resources) / median_time(larger_resources)
```

The reference here is one MPI rank with four threads, not the single-threaded serial executable. Report this baseline explicitly. Weak scaling does not compare the same physical domain at different resolutions.

## Resource layout and communication

At 16 total cores, compare explicit rank/thread combinations. For each resource layout, compare row and 2D decomposition, and blocking and overlap communication. Change one factor at a time when explaining a result.

Inspect local block dimensions and exposed boundary lengths. A layout with less communication is not guaranteed to be faster if it changes message overhead, local threading efficiency or memory placement.

## Planned figures

- Pigment field for a small demonstration case, with an optional velocity overlay.
- Maximum absolute serial/parallel difference for validation cases.
- Relative mass change for the validation cases.
- Analytical error versus resolution for the cosine case.
- Runtime, speedup and efficiency versus total cores for strong scaling.
- Runtime or weak-scaling efficiency versus resources for weak scaling.
- Fixed-core comparison of rank/thread layouts and communication strategies.

Use only retained, measured data. Label the machine, grid, steps, resource definition and timing scope. Separate observations from proposed explanations.
