# Model B: 3D flow with pressure projection and passive pigment

This is a self-contained CPU implementation of the Model B proposed in the
project's `docs/design_3d.md`. It solves **three-dimensional velocity, pressure
projection and pigment transport together at every time step**. It sits beside the 2D solvers, Model A and the standalone Poisson benchmark without changing them. It is an additional experiment, with separate result files.

## What is solved

In a closed rectangular basin, with constant density, the target equations are

\[
\partial_t\mathbf u+\nabla\cdot(\mathbf u\mathbf u)
=-\nabla(p/\rho)+\nu\Delta\mathbf u+\mathbf f,
\qquad \nabla\cdot\mathbf u=0,
\]

\[
\partial_t c+\nabla\cdot(\mathbf u c)=\kappa\Delta c,
\qquad \mathbf u=(u,v,w).
\]

All six walls are impermeable and free-slip: normal velocity is zero, tangential
velocity has zero normal derivative. Pigment has zero normal flux. There is no
buoyancy or feedback from pigment into momentum; the pigment is a passive scalar.
All axes and file indices increase in the positive coordinate direction, unlike
the downward row/layer ordering in Model A. No Model A field is read implicitly.

The default basin is `1 x 1 x 0.5`, with `32 x 32 x 16` cells. Pressure and pigment
are cell-centred, velocities are stored at east, north and top faces (MAC grid).
The default localized pigment patch is a Gaussian in normalized coordinates,
centred at `(0.3 Lx, 0.5 Ly, 0.5 Lz)` with widths `(0.08 Lx, 0.08 Ly, 0.08 Lz)`.
Flow starts from rest. Optional spin-up precedes pigment release.

The wind body force is

\[
f_x=\frac{F_0}{k_y}(1+\alpha\sin\omega t)\cos(k_y y)g(z),
\quad k_y=2\pi/L_y,\qquad f_y=f_z=0,
\]

with the surface profile

\[
g(z)=\frac{4\exp[-4(1-z/L_z)]}{1-\exp(-4)}.
\]

Its continuous depth average is one. Midpoint sampling need not have an exactly
unit discrete mean. `--profile uniform` sets `g=1` as a depth-independence control.
Surface forcing and the pressure correction produce depth-dependent horizontal
flow and nonzero vertical motion, despite the lack of an imposed vertical force.

## One step and the pressure boundary condition

1. Advance pigment with the current projected face velocity, conservative
   first-order upwind fluxes, central diffusion and forward Euler.
2. Evaluate staggered conservative central momentum advection, the seven-point
   viscous Laplacian and the wind force.
3. Predict velocity with Adams--Bashforth: Euler at the first flow step, AB2 at the
   second, AB3 thereafter. Force normal wall faces to remain zero.
4. Solve a **new pressure Poisson system on every step**:

\[
A\phi=b,\qquad A=-DG,\qquad b=-\frac{D\mathbf u^*}{\Delta t},
\qquad \partial_n\phi=0,\qquad \langle\phi\rangle=0.
\]

5. Correct the three velocity components:

\[
\mathbf u^{n+1}=\mathbf u^*-\Delta t\,G\phi.
\]

`D` and `G` are a compatible MAC divergence/gradient pair. The code removes the
global RHS mean to enforce the Neumann compatibility condition, reports that
removed mean, warm-starts from the preceding `phi`, and solves by matrix-free CG
in the mean-zero space. It verifies the **actual** residual `b-A phi` before
accepting convergence. The RHS has zero mean analytically because all normal
wall velocities are zero. A large removed mean would therefore indicate a bug,
not a physical source to be silently corrected.

The stopping test is `RMS(b-A phi) <= max(atol, tol*RMS(b))`, with defaults
`tol=1e-10`, `atol=1e-13`. The projection identity, including a nonzero removed
mean `m`, is `D u_new = -dt*(b-A phi + m)`. The default simulations check the
resulting divergence bound each step.

**Interpretation of pressure:** `phi` is the kinematic pressure potential that
enforces the projection over a time step. Under the AB predictor it is an
effective, time-weighted pressure contribution. This package does not claim a
third-order-accurate instantaneous physical pressure reconstruction. Multiply a
kinematic pressure by a specified density to express it in pressure units, but
this does not change that time-discretization qualification. The pigment scheme
is first order in time and upwind advection is first order in space. AB3 alone
does not make the complete model third order.

The old standalone benchmark has node-centred **Dirichlet** conditions. Its
measured results remain useful for that benchmark, but neither its old RHS nor
its boundary stencil can directly serve as this pressure solve. In this model,
the number of pressure degrees of freedom is `nx*ny*nz`, with one constant null
mode removed. The implied matrix has
`7*N - 2*(nx*ny + nx*nz + ny*nz)` nonzero entries; it is never explicitly stored.

## Parallel computation

- `model_b_serial`: a separate build with no MPI/OpenMP, one CPU thread.
- `model_b_mpi`: MPI domain decomposition with OpenMP loops within each block.
- `--decomp slabs`: split y; `pencils`: split x/y; `blocks`: split x/y/z.
  Optional `--px`, `--py`, `--pz` constrain the process grid. Automatic choices
  minimize a face-communication count over feasible factorizations.
- Pigment and CG exchange face ghosts. Velocity exchanges run through x, y and
  z sequentially, including ghosts, so diagonal momentum stencils receive edge
  data. Exchanges are **blocking** in this first version.
- CG matrix application requires a halo exchange. Global dot products use
  `MPI_Allreduce`; mean removal and actual-residual verification add reductions.
- Memory is initialized by parallel first-touch. MPI runs outside OpenMP regions
  using `MPI_THREAD_FUNNELED`.

**Implemented here:** serial CPU, MPI/OpenMP CPU, the full coupled time loop,
verification, CSVs and single-node PBS jobs. **Not implemented here:** GPU
pressure projection, communication overlap, restart, or validated multi-node
benchmark placement. The existing A40 measurements must not be labelled as
Model B GPU results. No HPC jobs are submitted by building or unpacking.

## Build and run locally

Requires C++17, CMake, Python 3.8+, and MPI/OpenMP for the parallel executable.
There are no Python package dependencies. On macOS, `build.sh` finds an installed
Homebrew libomp. On Vanda it uses the modules already used by this project.

```bash
cd model_b
bash build.sh
python3 verify.py --output validation \
  --launcher 'mpirun --host localhost --oversubscribe --bind-to none'
python3 run_model_b.py --suite smoke --output results/local-smoke \
  --launcher 'mpirun --host localhost --oversubscribe --bind-to none' \
  --layouts serial,4x2
```

The oversubscription option above is for small correctness tests, not reported
scaling experiments. To build/run only the serial version:

```bash
MODELB_SERIAL_ONLY=1 bash build.sh build-serial
python3 verify.py --serial-only --build-dir build-serial --output validation-serial
python3 run_model_b.py --build-dir build-serial --layouts serial --output results/serial
```

An example individual run is

```bash
./build/model_b_serial --nx 32 --ny 32 --nz 16 --steps 100 --spinup 0 \
  --nu 0.001 --kappa 0.0001 --force 0.8 --tol 1e-10 \
  --trace trace.csv --output final_field.csv
```

The time step is fixed over the entire run. Unless explicitly set with `--dt`,
it uses `--speed-bound 2` in each direction and conservative advective/diffusive
limits. The actual velocities are checked every step. A violation or CG
nonconvergence produces a nonzero exit code; the program does not silently
change dt (which would invalidate the fixed-step AB coefficients).

## Run on Vanda

From this folder in a checkout of the repository on Vanda:

```bash
cd model_b
qsub -P <project> vanda_smoke.pbs
```

`<project>` is your Vanda project (the scripts carry no `#PBS -P` line). The smoke job requests
one node, 16 CPU cores, 16 GB, 30 minutes. It first runs the verification suite,
then a 100-step serial/4-rank x 2-thread complete simulation on 32 x 32 x 16
cells. It returns nonzero if compilation, any check or any simulation fails.

After this job passes, a modest numerical experiment is available:

```bash
qsub -P <project> vanda_run.pbs
```

It requests one node, 32 CPU cores, 32 GB, one hour. It reruns the correctness
gate for its own build, then runs 64 x 64 x 32 cells, 100 steps, no spin-up,
4 MPI ranks x 8 threads, one warm-up and three measured repetitions. This is a
full-model cost measurement, not yet a scaling campaign. Both scripts put builds
and results in unique job-ID directories to avoid concurrent-job collisions.
`mpi_vanda.sh` follows the project's Open MPI 5 prefix/environment setup and
binds production threads to CPU cores. It deliberately rejects multi-node jobs.

Useful checks (replace the example job ID):

```bash
qstat -fx JOB_ID
cat results/JOB_ID/validation/status.json
cat results/JOB_ID/smoke/summary.csv
# For vanda_run.pbs instead:
cat results/JOB_ID/model_b/summary.csv
```

Copy a result to your own computer from a terminal there:

```bash
scp <user>@vanda.nus.edu.sg:<checkout>/model_b/results/JOB_ID/model_b/summary.csv .
```

## Outputs, timing and experimental interpretation

- `validation/validation.csv`: mathematical and parallel correctness checks.
- `raw.csv`: one row per run, phase/repetition, MPI/thread counts, time step,
  pressure iteration counts, divergence, mass, vertical velocity and timings.
- `summary.csv`: measured-run medians and min/max loop times. In benchmark mode,
  when a matching serial layout is present, also reports speedup and efficiency.
- `*_trace.csv` (smoke): per-step pressure residual, iteration count, divergence,
  concentration bounds, mass and maximum vertical velocity.
- `*_field.csv` (smoke): one row per cell, coordinates, pigment, `phi`, and the
  three high-face velocities. Velocities are at faces, not at the listed centre
  coordinates. Row order may differ by decomposition; compare by `(x,y,z)`.
- `.log`: exact launch command and program diagnostics; `setup.json` records the
  run settings. `status.json` is written PASS only after the suite completes.

`loop_s` is the maximum wall time across ranks of the complete coupled loop.
It includes pigment, momentum, pressure projection and per-step stability and
divergence checks. It excludes allocation/setup, spin-up, and file writes.
Smoke runs additionally collect per-step trace diagnostics in that timed loop;
benchmark runs disable them and full-field output. Do not compare their timings
without acknowledging that difference. `end_to_end_s` includes launch and writes.
`poisson_s` is a separately reduced maximum and is nested within the loop time;
`halo_s` includes halo packing/wall fills and `reduction_s` includes diagnostic
reductions. These component timings overlap and are not an additive breakdown.

Existing one-node/four-node/A40 measurements solve one fixed Dirichlet system.
New Model B measurements solve a changing Neumann pressure system once per time
step plus evolve the velocity and pigment. Use separate tables and definitions;
do not reuse an old serial reference to calculate the new model's speedup.

The included small-grid local validation is evidence for the implemented
operators and parallel consistency. It does not establish mesh independence,
long-time stability at all parameters, turbulent accuracy, instantaneous pressure
accuracy, or Vanda performance. Larger scientific studies need their own checks.

## Numerical references

- Project `docs/design_3d.md`, Model B proposal (local copy inspected on 16 Nov 2026).
- [AMReX-Hydro projection methods](https://amrex-fluids.github.io/amrex-hydro/docs_html/Projections.html)
  explains compatible discrete projection operators and pressure potentials.

The implementation is self-contained and does not depend on AMReX.

## Large runs: 1024^3 strong scaling and weak scaling to 256 cores (added 16 Nov 2026)

Four jobs, each with its own build and the 69-check gate, submitted from this folder
with `qsub -P <project> <job>`. All runs use the unit cube, cubic cells, one thread
per rank, `--tol 1e-8` (a relative 1e-10 is at the rounding limit of a 1024^3 solve)
and `--progress 100`, which prints the solve time every 100 CG iterations so that a
run stopped by the wall time still gives a time per iteration.

| Job | Nodes, wall time | Runs |
| --- | --- | --- |
| `big_strong_p1.pbs` | 1, 12 h | 1024^3, 1 rank, one time step |
| `big_strong_p2_p4.pbs` | 1, 12 h | 1024^3, 4 and 2 ranks, one time step |
| `big_strong_node.pbs` | 1, 8 h | 1024^3, 64 (3 repeats), 32, 16 and 8 ranks, one time step |
| `big_four_nodes.pbs` | 4, 8 h | serial against 256 ranks on 48^3 (fields compared); 1024^3 at 256 and 128 ranks (3 repeats); weak scaling, 2 steps, 1 warm-up and 2 repeats per point |

Weak scaling keeps about 2^22 cells per rank: 162^3, 204^3, 256^3, 322^3, 406^3,
512^3, 646^3, 812^3 and 1024^3 cells on 1, 2, 4, 8, 16, 32, 64, 128 and 256 ranks.
The CG iterations per step grow with the grid, so compare the time per CG iteration.

`mpi_vanda_big.sh` places the ranks as the 2D experiments do: up to 32 ranks by core,
more ranks 32 per socket (64 per node). Results are in `results/<job id>/<point>/`
(`raw.csv`, `summary.csv`, one log per run).
