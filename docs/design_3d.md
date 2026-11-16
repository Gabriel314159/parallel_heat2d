# Three-dimensional model

**Status, 16 November 2026.** The extension to three dimensions is done in two steps, and a standalone 3D Poisson benchmark sits beside them.

- **Model A, implemented:** the pigment in three dimensions, carried by the 2D wind-driven flow of `design.md`, the same in every layer and with no vertical velocity. Programs `heat3d_serial` and `heat3d_parallel`; see the next section.
- **Model B, proposal of 3 November, not implemented:** the 3D flow with a projection method, the 3D pressure Poisson equation and the GPU offload. It is described from "Why and what" to the end of this document, and its decisions are still open.
- **3D Poisson benchmark, implemented and measured:** `poisson3d/` solves $-\Delta\phi=f$ with the seven-point operator and Dirichlet walls by CG, in serial, with MPI + OpenMP (x and y split, as in model A's pencils) and on GPUs. It is not part of model A, whose flow still solves the 2D Poisson equation; it serves the 2D/3D Poisson comparison. See `poisson3d/README.md`.

## Model A: the 3D pigment in the 2D flow (implemented)

### What it is

The pigment fills the basin $[0,l_x]\times[0,l_y]\times[0,l_z]$ and obeys the 3D advection–diffusion equation

$$
\frac{\partial c}{\partial t} + \nabla\cdot(\mathbf u c) = \kappa\nabla^2 c,
\qquad
\mathbf u = \big(u(x,y,t),\ v(x,y,t),\ 0\big),
$$

where $u$ and $v$ are the 2D wind-driven flow of `design.md`, computed exactly as in `heat_serial` and `heat_parallel`. The flow is the same in every layer and has no vertical velocity. Its divergence is $\partial u/\partial x + \partial v/\partial y = 0$, as in 2D, so every property of the 2D pigment carries over: the total pigment is conserved, a constant field stays constant, and no new extrema appear. Vertically, the pigment only diffuses.

Model A gives the parallel parts of the 3D problem (a 3D field, a seven-point stencil, faces that are 2D slabs, decompositions of a 3D grid) with the flow code and its Poisson solver unchanged. Model B would add the 3D flow on top.

### Grid and storage

The grid is cell-centred, as in 2D, with a third index for the layer: $(k, i, j)$ = (layer, row, column). The layer $k$ counts down from the surface, as the row $i$ counts down from the north wall, so the centre of cell $(k,i,j)$ is at

$$
x = (j-\tfrac12)\Delta x, \qquad y = l_y - (i-\tfrac12)\Delta y, \qquad z = l_z - (k-\tfrac12)\Delta z,
\qquad \Delta z = l_z/n_z .
$$

`Field3D` (`src/common3d/field3d.hpp`) stores one block with a ghost layer on all six sides in one array, layer by layer and row by row, with $x$ fastest. So `row(k, i)[j]` walks along $x$ exactly as `Field::row(i)[j]` does in 2D. The face velocities stay the 2D `FaceVelocities`, read by every layer.

The defaults keep the 2D grid (256 × 256 cells on 1 × 1) and add $n_z = 32$ layers over $l_z = 0.125$, so the cells are cubes.

### Pigment update

The 2D update gains the fluxes through the top and bottom faces:

$$
c^{n+1}_{k,i,j} = c^n_{k,i,j}
- \frac{\Delta t}{\Delta x}\big(F_{k,i,j} - F_{k,i,j-1}\big)
- \frac{\Delta t}{\Delta y}\big(G_{k,i,j} - G_{k,i+1,j}\big)
- \frac{\Delta t}{\Delta z}\big(H^{\text{top}}_{k,i,j} - H^{\text{top}}_{k+1,i,j}\big).
$$

$F$ and $G$ are the fluxes of the 2D update (`face_flux` of `common/stencil.hpp`, through the east and north faces), and $H^{\text{top}}_{k,i,j}$ is the flux through the top face of cell $(k,i,j)$, positive upwards. It is again `face_flux`, with the cell below the face on its minus side and $w = 0$:

$$
H^{\text{top}}_{k,i,j} = 0\cdot c_{k-1,i,j} - \frac{\kappa}{\Delta z}\big(c_{k-1,i,j} - c_{k,i,j}\big).
$$

`update_cell_3d` (`src/common3d/stencil3d.hpp`) takes $w$ as an argument, so a 3D flow could use it; model A passes 0. Both cells that share a face compute its flux with the same arguments in the same order, so the total pigment is conserved to rounding, as in 2D.

**Walls.** All six walls are no-flux. The ghosts across them mirror the interior: the four side walls of a block, layer by layer, as `fill_wall_ghosts` does in 2D, and the surface and the bottom, which every block touches (`fill_wall_ghosts_3d` in `src/common3d/physics3d.cpp`).

**Reduction to 2D.** When the layers above and below hold the same value as the cell (one layer between two mirrored ghosts, or a field that does not change with depth), both vertical fluxes are exactly zero, and the new value equals the 2D update exactly. So a 3D run with one layer, or with the 2D drop in every layer (`--blob-z column`), must reproduce the 2D pigment exactly in every layer; `tests/test_3d.cpp` checks this.

### Initial conditions

- `pigment_drop`: a 3D Gaussian, $c_0 + A_0\exp\!\big(-|\mathbf r - \mathbf r_0|^2/(2\sigma^2)\big)$, centred at $(x_0, y_0, z_0)$ with $z_0$ = `--blob-z` (default $l_z/2$). With `--blob-z column`, the 2D drop in every layer.
- `constant`: $c_0$.
- `cosine_decay`: $A_0\cos(\pi x/l_x)\cos(\pi y/l_y)\cos(\pi z/l_z)$, without wind. It is an exact eigenvector of the scheme, with
  $$
  \lambda = 4\kappa\left(\frac{\sin^2(\pi\Delta x/2l_x)}{\Delta x^2} + \frac{\sin^2(\pi\Delta y/2l_y)}{\Delta y^2} + \frac{\sin^2(\pi\Delta z/2l_z)}{\Delta z^2}\right),
  $$
  so after $n$ steps the field is $(1-\lambda\Delta t)^n$ times the initial one.

### Time step

The rules of `design.md` with the vertical diffusion added to the pigment limit:

$$
a = \frac{u_{\max}}{\Delta x} + \frac{v_{\max}}{\Delta y},
\qquad
d = 2\kappa\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2}\right),
\qquad
n = \nu\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right),
$$

with $\Delta t\,(a+d) \le 1$, $\Delta t\,a \le 0.5$ and $\Delta t\,n \le 0.1$. There is no vertical velocity, and the flow, hence its viscous limit, stays 2D. Only the pigment limit is tighter than in 2D, so the 3D step is never larger than the 2D step for the same $x$-$y$ grid.

### The flow, live or frozen

Each step advances the 2D flow exactly as in 2D: the same code (`SerialFlow` in `src/serial/serial_flow.hpp`, `ParallelFlow` in `src/parallel/parallel_flow.hpp`, moved there from the 2D solvers unchanged), with one 2D Poisson solve.

With `--freeze-flow`, the flow is advanced only during the spin-up, and the pigment then moves in that fixed flow. No Poisson solve runs in the timed loop, which then measures the 3D pigment alone. This is what the 3D experiments use: with the flow live, the 2D CG, already measured in experiments 1 to 5, takes most of the time again.

### Decomposition and halo exchange

Every rank owns whole columns of water: the ranks split $x$ and $y$ and never $z$.

- `--decomp slabs` splits $y$ only, like the 2D rows;
- `--decomp pencils` splits $x$ and $y$, like the 2D blocks (`--py`, `--px`, or `MPI_Dims_create`).

The 2D layout code (`make_layout`) is used as it is. A rank therefore has the same four neighbours as in 2D, and its 2D flow block and halo exchanges are those of `heat_parallel`. Splitting $z$ as well (`blocks`) is left to model B: ranks with the same columns of water would all need the same 2D flow, which they would have to compute several times over or receive from one another.

Each pigment face a rank exchanges is a 2D slab: one row or one column in every layer (`Halo3D`, `src/parallel3d/halo3d.cpp`).

- The north and south faces (row 1 or $n_y$ of every layer) are $n_z$ blocks of $n_x$ contiguous values, one layer apart. They are sent straight from the field with one MPI vector datatype.
- The west and east faces (column 1 or $n_x$ of every row of every layer) are packed into buffers of $n_z \times n_y$ values, as the 2D columns are.

The seven-point stencil reads no edge or corner ghosts, so faces are enough. With `--exchange overlap`, the messages are started, the wall ghosts are filled, the strict interior of the 2D rule is updated in every layer, the messages are waited for, and the remaining cells of every layer are updated exactly once.

For the 3D benchmark problem $P_3$ (1024 × 1024 × 64 cells, see "Experiments" below), the values a rank sends per step against the cells it owns are:

| Ranks | Decomposition | Process grid ($y \times x$) | Cells per rank ($x \times y \times z$) | Halo values / owned values |
| --- | --- | --- | --- | --- |
| 64 | slabs | 64 × 1 | 1024 × 16 × 64 | 0.125 |
| 64 | pencils | 8 × 8 | 128 × 128 × 64 | 0.031 |
| 256 | slabs | 256 × 1 | 1024 × 4 × 64 | 0.5 |
| 256 | pencils | 16 × 16 | 64 × 64 × 64 | 0.063 |

These are upper bounds that ignore the walls.

### Output

- `--output PATH` writes the whole 3D field: CSV with $n_z n_y$ rows of $n_x$ values (layer by layer, surface first; within a layer, the north row first, as in 2D), or `.npy` of shape $(n_z, n_y, n_x)$. The parallel solver gathers the field only when it is written.
- Snapshots keep the 2D folder of `design.md`, with the surface layer as the pigment (`field_NNNNNN.npy`), so the 2D plotting scripts work unchanged, and add a vertical section through the row of the blob centre (`section_NNNNNN.npy`, $n_z \times n_x$) and `meta3d.json`.
- The `RESULT` line is the 2D one with the 3D fields added: `model`, `nz`, `lz`, `dz`, `blob_z`, `flow` (`live` or `frozen`), `pigment_s` (the pigment updates of the timed loop, halos included) and `pigment_halo_s` (their halo exchanges), each the largest over the ranks. `cell_updates_per_s` counts all $n_x n_y n_z$ cells.

### Correctness checks

- `tests/test_3d.cpp`: the 3D storage, options and time step; **one layer equals the 2D pigment exactly**, and so does **every layer of a depth-independent field**; a 3D blob in the flow (mass conserved to rounding, no new extrema, symmetric about mid-depth); a constant field under the flow; the frozen flow; the 3D cosine decay against the exact discrete solution and its convergence to the continuous one (at least second order; the error falls about 16-fold per refinement here, since the time and space errors partly cancel at this step); the output layouts.
- `tests/test_halo3d.cpp` (MPI, 1, 4 and 6 ranks): after a blocking or a nonblocking exchange, every face ghost next to another rank holds its neighbour's value in every layer, for slabs and pencils, and no wall ghost is touched.
- `experiments/correctness_3d.json` (experiment 0 for model A, run by `scripts/check_correctness.py`): 14 cases with serial checks and serial/parallel comparisons of the whole 3D field, with slabs and pencils, blocking and overlapping exchange, up to 6 ranks × 2 threads, uneven grids, one layer, 3 × 3 × 2 cells on 4 ranks, the periodic wind, the frozen flow, the depth-independent drop, the constant field and the cosine decay.

### Experiments

$P_3$ (`experiments/benchmark_problem_3d.json`) is the 2D problem $P$ with 64 layers over a depth of 0.0625 (cubic cells, 67 million of them), the blob at mid-depth and the flow frozen from rest. Three experiments use it, with the runner of the 2D experiments (`"model": "3d"` in the file selects the 3D programs, and the raw CSV files get the columns of `results/raw/metrics_columns_3d.csv`):

| File | Question | What changes |
| --- | --- | --- |
| `strong_scaling_3d.json` | How much faster do more ranks make $P_3$? | 1 to 256 ranks, one thread each, pencils, blocking |
| `communication_3d.json` | Do slabs or pencils, blocking or overlap matter when every face is 64 layers deep? | slabs or pencils × blocking or overlap, 256 ranks |
| `weak_scaling_3d.json` | Does the time stay the same with 256 × 256 × 64 cells per rank? | 1, 4, 16, 64 and 256 ranks |

The jobs are `hpc/jobs/exp3d_one_node.pbs` and `hpc/jobs/exp3d_four_nodes.pbs`. As for $P$, time one serial run on Vanda first and adjust $n_z$ or the number of steps in the one problem file.

### Code

| Path | Contents |
| --- | --- |
| `src/common3d/` | `Field3D`, the 3D stencil, the 3D options and time step (`Config3D`, on top of the 2D `Config`), initial conditions, walls, mass and cosine reference, and the 3D output. No MPI or OpenMP. |
| `src/serial3d/` | `heat3d_serial`. |
| `src/parallel3d/` | `heat3d_parallel`: the 3D halo exchange, the gathers for output and the solver. |
| `src/serial/serial_flow.hpp`, `src/parallel/parallel_flow.hpp`, `src/parallel/parallel_tools.hpp` | The 2D flow of each solver and the MPI helpers, moved out of the 2D solvers unchanged so that both models share them. |

The 2D programs are unchanged in behaviour; their tests and correctness cases pass as before.

## Why and what

*This section and all after it describe model B, the proposal of 3 November.*

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
