# Validation and performance plan

This document defines the experiment protocol. Completed local numerical and pipeline checks are recorded in [verification.md](verification.md). Atlas/Vanda performance measurements remain pending.

## Correctness

### Serial and parallel equivalence

For each small case, compare the final full matrix from `heat_serial` with each selected parallel strategy. Use identical global grid dimensions, geometry, diffusivity, boundary/initial conditions, time step, step count and precision. Compare with an absolute and relative tolerance, initially `1e-10` each, and record the observed error.

Test one rank, several ranks, multiple threads, non-divisible dimensions and very small local blocks. Include both blocking and overlap modes. Check the physical boundary values separately.

### Constant field

Set all interior and boundary values to the same constant. Every numerical update should retain that constant up to floating-point roundoff. This checks initialization, boundary handling, halo updates and the update formula.

### Analytical sine decay

For zero-temperature boundaries on `[0,lx] x [0,ly]`, initialize:

```text
T(x,y,0) = A * sin(pi*x/lx) * sin(pi*y/ly)

T(x,y,t) = A * sin(pi*x/lx) * sin(pi*y/ly)
           * exp(-alpha*pi^2*(1/lx^2 + 1/ly^2)*t)
```

Evaluate at the actual `t = steps*dt`. A finite-difference solution has discretization error; it is not expected to match this continuous solution to the serial/parallel equivalence tolerance. Run a separate mesh-refinement study at a common physical final time, with stable time steps, and check that errors decrease. Record both time and spatial resolution.

The correctness manifest includes one sine case for baseline equivalence and analytical error reporting. It does not by itself define the complete refinement study; add the common-final-time resolution series when implementing that study.

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

Keep a 512-by-512 interior block and four threads per rank, then use 1, 4 and 16 ranks with square process grids. Keep `dx` and `dy` fixed and grow the physical domain according to the configuration rules. Thus the automatic `dt`, step count and work per core remain comparable.

```text
weak-scaling efficiency = median_time(reference_resources) / median_time(larger_resources)
```

The reference here is one MPI rank with four threads, not the single-threaded serial executable. Report this baseline explicitly. Weak scaling does not compare the same physical plate at different resolutions.

## Resource layout and communication

At 16 total cores, compare explicit rank/thread combinations. For each resource layout, compare row and 2D decomposition, and blocking and overlap communication. Change one factor at a time when explaining a result.

Inspect local block dimensions and exposed boundary lengths. A layout with less communication is not guaranteed to be faster if it changes message overhead, local threading efficiency or memory placement.

## Planned figures

- Temperature field for a small demonstration case.
- Maximum absolute serial/parallel difference for validation cases.
- Analytical error versus resolution for the sine case.
- Runtime, speedup and efficiency versus total cores for strong scaling.
- Runtime or weak-scaling efficiency versus resources for weak scaling.
- Fixed-core comparison of rank/thread layouts and communication strategies.

Use only retained, measured data. Label the machine, grid, steps, resource definition and timing scope. Separate observations from proposed explanations.
