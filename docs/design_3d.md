# Three-dimensional model: design proposal

**Status:** proposal for review, 3 November 2026. Nothing here is implemented. The decisions still open are listed at the end; until they are settled, the 2D model of `design.md` stays the model we measure.

## Why and what

We extend the whole model to three dimensions: the pigment and the wind-driven flow. The Poisson solve, which takes over 99% of the run time in 2D, becomes a three-dimensional pressure Poisson equation, and in experiment 3 it also runs on GPUs through OpenMP target offload.

The streamfunction of the 2D model has no simple counterpart in 3D, so the flow changes formulation: from vorticity and streamfunction to velocity and pressure, advanced with a projection method. Everything else carries over.

| | 2D model (`design.md`) | 3D model |
| --- | --- | --- |
| Pigment | finite volumes, upwind advection, five-point diffusion, forward Euler | the same, with a third direction (seven-point diffusion) |
| Flow unknowns | $\zeta$, $\psi$ at cell corners | $u, v, w$ on cell faces, pressure $\phi$ at cell centres |
| Advection of the flow | Arakawa Jacobian | central differences in divergence form |
| Time stepping of the flow | AB3 | AB3, followed by a projection |
| Poisson equation | $-\nabla^2\psi = \zeta$, $\psi = 0$ on the walls, five points | $-\nabla^2\phi = -\nabla\cdot\mathbf u^*/\Delta t$, $\partial\phi/\partial n = 0$ on the walls, seven points |
| Poisson solver | CG, matrix-free, warm start | the same, with the mean removed from the right-hand side |
| Divergence of $\mathbf u$ | zero to rounding | zero to the CG tolerance |
| Decompositions | rows, 2D blocks | slabs, pencils, 3D blocks |

## Model

The basin is $[0,l_x]\times[0,l_y]\times[0,l_z]$, with $z$ pointing up and the surface at $z = l_z$. The pigment obeys

$$
\frac{\partial c}{\partial t} + \nabla\cdot(\mathbf u c) = \kappa\nabla^2 c,
\qquad
\mathbf u = (u, v, w).
$$

The flow obeys the incompressible Navier–Stokes equations,

$$
\frac{\partial \mathbf u}{\partial t} + \nabla\cdot(\mathbf u\mathbf u) = -\nabla\phi + \nu\nabla^2\mathbf u + \mathbf f,
\qquad
\nabla\cdot\mathbf u = 0,
$$

where $\phi$ is the pressure divided by the density. It has no equation of its own: it is whatever keeps $\mathbf u$ divergence-free.

### Wind forcing

The wind pushes the water along $x$, more strongly near the surface:

$$
\mathbf f = \big(f_x(y,z,t),\,0,\,0\big),
\qquad
f_x = \frac{F_0}{k}\,\big(1 + \alpha\sin(\omega t)\big)\cos(ky)\,g(z),
\qquad
k = \frac{2\pi}{l_y}.
$$

Its curl about the vertical is

$$
\frac{\partial f_y}{\partial x} - \frac{\partial f_x}{\partial y} = F_0\,\big(1 + \alpha\sin(\omega t)\big)\sin(ky)\,g(z),
$$

which is the 2D forcing $F(y,t)$ times $g(z)$. So every layer is driven towards the double gyre of the 2D model.

The depth profile $g$ has mean 1 over the depth, so that the depth-averaged forcing equals the 2D forcing:

- `surface` (default): $g(z) = \dfrac{l_z}{h}\,\dfrac{e^{-(l_z - z)/h}}{1 - e^{-l_z/h}}$, which decays over a depth $h$ (default $h = l_z/4$);
- `uniform`: $g = 1$. A flow that starts independent of $z$ then stays independent of $z$, with $w = 0$, which gives a test (see "Correctness checks").

### Walls

All six walls are free-slip and no-flux:

- no water crosses a wall: the normal velocity is 0;
- the wall exerts no drag: the tangential velocity has zero normal derivative;
- no pigment crosses a wall, as in 2D.

These are the 3D form of the 2D walls ($\psi = \zeta = 0$). The parameters $U$, $\mathrm{Re}$, $\nu$ and $F_0$ are chosen as in 2D.

## Grid

The grid is cell-centred, as in 2D, with a third index for the layer:

- $(k, i, j)$ = (layer, row, column). The layer $k$ counts down from the surface, as the row $i$ counts down from the north wall.
- The cell centre is at $x = (j-\tfrac12)\Delta x$, $y = l_y - (i-\tfrac12)\Delta y$ and $z = l_z - (k-\tfrac12)\Delta z$.
- $c$ and $\phi$ are stored at the cell centres, with one ghost layer on every side.

The velocity lives on the faces (a staggered, or MAC, grid). Each component is stored in a field of the same shape as $c$. Entry $(k,i,j)$ holds the face on the high-index side of the cell:

$$
u_{k,i,j}\ \text{on the east face}\ x = j\Delta x,
\qquad
v_{k,i,j}\ \text{on the south face}\ y = l_y - i\Delta y,
\qquad
w_{k,i,j}\ \text{on the bottom face}\ z = l_z - k\Delta z.
$$

This is the rule of the 2D corner fields, which hold the south-east corner of each cell. The wall faces of index $n_d$ lie in the last owned layer of the blocks on that wall. The wall faces of index 0 are in the ghost layer. All wall faces hold 0.

The pigment needs no interpolation: its upwind fluxes use these face velocities directly, as the 2D pigment uses $U_e$ and $V_n$.

## Discretization

### Pigment

The pigment update is the 2D update with a third pair of faces. With the total flux through each face (upwind advection plus central diffusion), it is

$$
c^{n+1}_{k,i,j} = c^n_{k,i,j}
- \frac{\Delta t}{\Delta x}\big(F_{k,i,j} - F_{k,i,j-1}\big)
- \frac{\Delta t}{\Delta y}\big(G_{k,i-1,j} - G_{k,i,j}\big)
- \frac{\Delta t}{\Delta z}\big(H_{k-1,i,j} - H_{k,i,j}\big),
$$

where $F$, $G$ and $H$ are the fluxes through the east, south and bottom faces, positive along $x$, $y$ and $z$. The no-flux ghost fill and forward Euler are unchanged, and the total pigment is still conserved to rounding.

### Momentum

We write the right-hand side of the momentum equation without the pressure:

$$
\mathbf R = -\nabla\cdot(\mathbf u\mathbf u) + \nu\nabla^2\mathbf u + \mathbf f .
$$

Each component is evaluated at its own faces. We use two operators along each direction $d$:

- $\overline{a}^{d}$, the average of two neighbours a half cell apart;
- $\delta_d a$, their difference divided by the cell width.

The advection term is then the divergence form on the staggered grid. For $u$ it reads

$$
\big(\nabla\cdot(\mathbf u u)\big)_u
= \delta_x\big(\overline{u}^{x}\,\overline{u}^{x}\big)
+ \delta_y\big(\overline{v}^{x}\,\overline{u}^{y}\big)
+ \delta_z\big(\overline{w}^{x}\,\overline{u}^{z}\big),
$$

and similarly for $v$ and $w$. It is second order. For a velocity that is discretely divergence-free, it conserves momentum and kinetic energy (Morinishi et al., *J. Comput. Phys.* 143, 1998). This is the role the Arakawa Jacobian plays in 2D.

The viscosity is the seven-point Laplacian of each component. At free-slip walls, the tangential components are mirrored into the ghost layer. The forcing is evaluated at the $u$ faces, at the old time.

The stencils read diagonal neighbours: $u$ at an east face needs $v$ on the south faces of the cells on both sides of it, one of them in the neighbouring column. So the velocity halos must include the edge ghosts. The 2D corner fields have the same need (see "Halo exchange").

### Time stepping and projection

One flow step has three parts.

1. **Predict.** Advance without the pressure, with AB3 as in 2D (Euler at the first step, AB2 at the second):

   $$
   \mathbf u^* = \mathbf u^n + \Delta t\left(\tfrac{23}{12}\mathbf R^n - \tfrac{16}{12}\mathbf R^{n-1} + \tfrac{5}{12}\mathbf R^{n-2}\right).
   $$

2. **Solve for the pressure.** Find $\phi$ such that the corrected velocity is divergence-free:

   $$
   D\,G\,\phi = \frac{D\,\mathbf u^*}{\Delta t}.
   $$

   $D$ is the divergence of a cell, from its six faces, and $G$ is the gradient on a face, from its two cells. $G$ is 0 on the wall faces.

3. **Correct.** $\mathbf u^{n+1} = \mathbf u^* - \Delta t\,G\phi$.

Every term is explicit, so the projection adds no splitting error, and the scheme is third order in time like the 2D vorticity. There is one Poisson solve per step, as in 2D.

### Poisson equation for φ

$D\,G$ is the seven-point Laplacian, with the neighbour across a wall left out. In a cell next to the west wall, for example, the $x$ part is $(\phi_E - \phi)/\Delta x^2$. We solve

$$
A\phi = b,
\qquad
A = -D\,G,
\qquad
b = -\frac{D\,\mathbf u^*}{\Delta t},
$$

with the CG of `design.md`, matrix-free and warm-started from the previous step's $\phi$.

Two things differ from the 2D Poisson equation.

- **$A$ is only semi-definite.** Adding a constant to $\phi$ changes nothing, so $A$ has the constants as its null space. CG still converges when $b$ has zero sum, which holds exactly in exact arithmetic: the sum of all cell divergences is the net flow through the walls, which is 0. In floating point the sum is only close to 0. So we subtract the mean of $b$ before each solve. This costs one more global sum per solve, not per iteration.
- **The divergence is zero only to the tolerance.** With the residual $r = b - A\phi$, the corrected velocity has divergence $D\,\mathbf u^{n+1} = -\Delta t\, r$. The stopping test $\lVert r\rVert_2 \le \text{tol}\,\lVert b\rVert_2$ therefore gives

  $$
  \lVert D\,\mathbf u^{n+1}\rVert_2 \le \text{tol}\,\lVert D\,\mathbf u^*\rVert_2 .
  $$

  The total pigment is still conserved to rounding, because the flux form telescopes whatever the velocity. A constant pigment, however, now stays constant only to about $\text{tol}$, not to rounding.

Each CG iteration still needs one halo exchange and two global sums. $A$ reads only the six face neighbours, so this exchange needs no edge ghosts.

### One time step

1. Update the pigment from $c^n$ to $c^{n+1}$ with the face velocities $\mathbf u^n$.
2. Compute $\mathbf R^n$ on all faces and predict $\mathbf u^*$.
3. Compute $b$ from $D\,\mathbf u^*$, remove its mean, and solve for $\phi$, starting from the previous $\phi$.
4. Correct to $\mathbf u^{n+1}$.

## Time step

The rules of `design.md` gain a third term:

$$
a = \frac{u_{\max}}{\Delta x} + \frac{v_{\max}}{\Delta y} + \frac{w_{\max}}{\Delta z},
\qquad
d = 2\kappa\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}\right).
$$

The limits are:

- pigment: $\Delta t\,(a + d) \le 1$;
- centred flow advection with AB3: $\Delta t\,a \le 0.5$;
- viscosity: $\nu\,\Delta t\,(\Delta x^{-2} + \Delta y^{-2} + \Delta z^{-2}) \le 0.1$.

The most negative eigenvalue of the viscous term is now $-4\nu\sum_d \Delta d^{-2}$, so the AB3 bound of $6/11$ again gives about $0.136$.

## Correctness checks

The 2D checks carry over, and the projection adds its own.

1. **Cosine decay.** With no wind, $A_0\cos(\pi x/l_x)\cos(\pi y/l_y)\cos(\pi z/l_z)$ is an exact eigenvector of the discrete scheme. $\lambda$ gains a third term $\tfrac{4\kappa}{\Delta z^2}\sin^2\!\big(\tfrac{\pi\Delta z}{2l_z}\big)$, and the solution after $n$ steps is $(1-\lambda\Delta t)^n$ times the initial field.
2. **Poisson solver.** The same cosine modes are eigenvectors of $A$ with Neumann walls. A right-hand side equal to one mode must give that mode divided by its eigenvalue, up to the tolerance. A constant right-hand side must be removed by the mean subtraction.
3. **Projection.** After every step, $\max\lvert D\,\mathbf u\rvert$ stays at the level of $\text{tol}\,\lVert D\,\mathbf u^*\rVert$; the run reports it at every check interval.
4. **Conservation.** The total pigment is constant to rounding.
5. **Depth independence.** With `uniform` forcing and a flow starting at rest, all layers stay equal and $w = 0$, to rounding. This must also hold when the decomposition splits $z$.
6. **Serial against parallel** (experiment 0). Every decomposition, exchange mode and thread count gives the serial answer, to rounding and the CG tolerance.
7. **GPU against CPU.** The same run with the solve offloaded gives the same answer, to the same level.

## Domain decomposition

There are three decompositions (`--decomposition`):

- `slabs` split $y$ only, like the 2D rows;
- `pencils` split $x$ and $y$, like the 2D blocks, and keep whole columns of water in one rank;
- `blocks` split all three directions.

Each direction is split as in 2D, with the remainder going to the first positions. The process grid is the factorisation of $p$ that sends the fewest halo values, found by trying every factorisation. `MPI_Dims_create` is not used, because it ignores the grid's shape. The Cartesian communicator is created without reordering, as in 2D.

For the grid proposed below, $128 \times 256 \times 64$ cells in $x \times y \times z$, the counts are:

| Ranks | Decomposition | Process grid ($x \times y \times z$) | Cells per rank | Halo values / owned values |
| --- | --- | --- | --- | --- |
| 64 | slabs | 1 × 64 × 1 | 128 × 4 × 64 | 0.50 |
| 64 | pencils | 4 × 16 × 1 | 32 × 16 × 64 | 0.19 |
| 64 | blocks | 4 × 8 × 2 | 32 × 32 × 32 | 0.19 |
| 256 | slabs | 1 × 256 × 1 | 128 × 1 × 64 | 2.00 |
| 256 | pencils | 8 × 32 × 1 | 16 × 8 × 64 | 0.38 |
| 256 | blocks | 4 × 16 × 4 | 32 × 16 × 16 | 0.31 |

The halo counts are upper bounds that ignore the walls. On 256 ranks the slabs are one cell thick: every cell is on an edge, and overlap has nothing to hide behind.

### Halo exchange

A block has up to six neighbours.

- **Pigment and CG.** They read only the six face neighbours. All six exchanges can be in flight at once, so they can be blocking or overlapping (`--exchange`). For overlap, the block splits into an inner box and six face strips, each cell exactly once, as in `overlap_split`.
- **Velocity.** Its stencils also read edge ghosts. We exchange one direction after another ($x$, then $y$, then $z$), each message including the ghost layers that the previous directions have just filled. This is the two-phase corner exchange of 2D with a third phase: six messages instead of one to each of the 26 neighbours. It is blocking.

All faces are packed into contiguous buffers, also those that are contiguous in memory. This keeps one code path, and on a GPU it lets a whole face be copied in one transfer.

## GPU offload

We offload the CG solve, the part that dominates the run time. The rest of the step stays on the CPU threads of the rank. If the momentum update turns out to take a large share in 3D, it can be offloaded next in the same way.

- **One GPU per rank.** Each rank takes the device `local_rank mod n_devices`, with the local rank from `MPI_COMM_TYPE_SHARED`. A Vanda GPU node has two A40s, so the runs are 1 rank with 1 GPU and 2 ranks with 2 GPUs (`experiment_plan.md`, experiment 3).
- **Data stays on the GPU.** $r$, $p$ and $Ap$ are allocated on the device once (`omp target enter data`). Each solve copies $b$ and the initial $\phi$ to the GPU and $\phi$ back: two fields per step, against hundreds of iterations on the device.
- **Loops.** The three CG loops become `omp target teams distribute parallel for collapse(3)`, with `reduction(+ : ...)` for the dot products. The two global sums stay on the host: each needs one number from the device.
- **Halos.** A rank packs its faces on the GPU, copies the buffers to the host, exchanges them with MPI and copies the received buffers back. With one rank, every side is a wall and nothing is copied during the iterations. With two ranks, one face goes each way per iteration.
- **Build.** A CMake option (`HEAT_OFFLOAD`) compiles the same pragmas for the GPU with `nvc++ -mp=gpu` and the MPI of `nvhpc-hpcx-cuda12/26.1`. Without the option, a macro turns them into the CPU `omp parallel for` loops, and GCC builds the code as before. Both compilers ran an offloaded test on Vanda's A40s (`hpc/README.md`).

The A40 computes in double precision at 1/64 of its single-precision rate. CG, however, moves about 90 bytes per cell and iteration and does about 20 floating-point operations, so its memory bandwidth (about 700 GB/s) should limit it first. A CPU node has comparable bandwidth over its two sockets, so one GPU against one CPU node is the fair comparison, together with the 256-core layouts.

## Benchmark problem and experiments

**Proposed problem $P$:** a basin of $1 \times 2 \times 0.5$ with $128 \times 256 \times 64$ cells, so $\Delta x = \Delta y = \Delta z = 1/128$; 2.1 million cells, four times the 2D $P$; 100 steps from rest. With about 19 fields of doubles, it needs about 320 MB, which fits one rank, one node and one GPU.

Its cost is still an estimate. CG needs about $\sqrt{\kappa(A)}$ iterations, with $\kappa(A)$ the ratio of the largest to the smallest nonzero eigenvalue of $A$. In 2D, $256 \times 512$ cells take about 350 iterations, 1.7 times $\sqrt{\kappa}$. Here the smallest nonzero eigenvalue is that of the mode $\cos(\pi y/l_y)$, $\pi^2/l_y^2$, and the largest is about $12/\Delta x^2$, so $\sqrt{\kappa} \approx 280$. With the same factor, that is about 500 iterations per step. The 2D CG takes about 5 ns per cell and iteration on the laptop, which gives about 5 s per step and 8 to 10 minutes per serial run. We time one serial run first and adjust $n_z$ or the number of steps, as for the 2D $P$.

The experiments keep their questions:

| # | Change for 3D |
| --- | --- |
| 0 | add checks 1 to 7 above |
| 1, 2 | the same ranks and threads, on the 3D $P$ |
| 3 | the same seven CPU layouts on 256 cores (`blocks`), plus 1 GPU and 2 GPUs |
| 4 | `slabs`, `pencils` and `blocks` × blocking and overlap: six variants on 256 ranks |
| 5 | a fixed block per rank (to be set together with $P$), with the domain growing in $x$ and $y$ |

## Code

The 3D model gets its own executables, `heat3d_serial` and `heat3d_parallel`, next to the 2D ones, which stay unchanged as the fallback.

- **New:** a 3D field type with layer, row and column strides; the 3D stencils for pigment and momentum; the projection; a 3D CG on the same `PoissonComm` interface; the 3D decomposition and halo exchange; the offload macro and device assignment.
- **Shared with 2D:** configuration parsing (extended by $l_z$, $n_z$ and the depth profile), timers, the run report and `RESULT` line, the binding check, node counting and MPI error handling.
- **Output:** the surface layer and one vertical section through the middle of the basin, written as 2D field files, so that the existing plots work. A full 3D snapshot is optional.

## Schedule

The report is due on 13 November, which leaves ten days. The benchmark runner is still missing and must be written in the same period.

| Dates | Work |
| --- | --- |
| 4–5 Nov | serial 3D model with checks 1 to 5 |
| 6–7 Nov | MPI decompositions and OpenMP; experiment 0 on the laptop and on Vanda; benchmark runner |
| **7 Nov** | **Decision point:** if serial and parallel do not agree by the evening, we measure the 2D model and report the 3D model as partial work |
| 8 Nov | CG offload and a GPU test job |
| 9–11 Nov | the runs on Vanda (experiments 0 to 5) |
| 11–13 Nov | report |

## Decisions for review

1. **Walls:** free-slip on all six walls (proposed), or no-slip, which needs a different ghost rule for the tangential velocity.
2. **Forcing:** the `surface` profile with $h = l_z/4$ (proposed), or uniform over the depth.
3. **Problem $P$:** $128 \times 256 \times 64$ cells in a $1 \times 2 \times 0.5$ basin (proposed), or $n_z = 32$ for runs half as long. With $n_z = 32$, pencils and blocks send the same amount on 256 ranks, so experiment 4 would separate them less.
4. **Experiment 4:** keep the slabs, which are one cell thick on 256 ranks, or compare only pencils and blocks.
5. **The 2D model in the report:** keep it as the starting point and as a 2D baseline, or describe the 3D model alone.
6. **The source of the requirement:** the wording that asks for the 3D Poisson solver and the GPU run, to cite in the report's introduction.
