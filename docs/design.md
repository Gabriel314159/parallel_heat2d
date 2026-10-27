# Design notes

The exact numerical and parallel implementations are not discussed here. Please refer to `experiment_plan.md` for details.

## Model

We aim to solve a two-dimensional advection–diffusion model of a passive scalar $c(x,y,t)$, for example pigment concentration, on the designated rectangular bound $[0,l_x] \times [0,l_y]$:

$$
\frac{\partial c}{\partial t}
+ \nabla \cdot (\mathbf{u}c)
= \kappa \left(
\frac{\partial^2 c}{\partial x^2}
+ \frac{\partial^2 c}{\partial y^2}
\right)
$$

- $\kappa > 0$ is a constant diffusivity left to our choice and is prescribed in `defaults.json`.
- $\mathbf{u} = (u, v)$ is a **prescribed** velocity field. It is given as a function of position and time. The scalar does not affect the flow. The velocity is derived from a streamfunction $\psi(x,y,t)$:

  $$
  u = \frac{\partial \psi}{\partial y},
  \qquad
  v = -\frac{\partial \psi}{\partial x}.
  $$

  It is therefore exactly divergence-free, or incompressible in the sense of fluid dynamics, because the mixed derivatives cancel:

  $$
  \nabla\cdot\mathbf{u}
  = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}
  = \frac{\partial^2 \psi}{\partial x\,\partial y} - \frac{\partial^2 \psi}{\partial y\,\partial x}
  = 0.
  $$

- The EXACT flow of the model is the time-periodic double gyre, a standard test flow for chaotic mixing. With $X = 2x/l_x$, $Y = y/l_y$, we define

  $$
  \psi(x,y,t) = A\,\sin\!\big(\pi f(X,t)\big)\,\sin(\pi Y),\\
$$
$$
  f(X,t) = X+\varepsilon\sin(\omega t)\,(X^2 -2X)
  $$

  The amplitude $A$ sets the flow speed, $\varepsilon$ the oscillation strength, and $\varepsilon = 0$ gives the steady double gyre, namely $f(X,t)=X$.
  Here, two counter-rotating gyres fill the box, and for $\varepsilon>0$ the dividing line between them oscillates left and right with period $2\pi/\omega$, which stretches and folds the scalar.  Since $f(0,t)=0$ and $f(2,t)=2$ for all $t$, $\psi = 0$ on all four walls, so the box is a closed container. 
- Because $\nabla\cdot\mathbf{u} = 0$, the advection term can be written as $\nabla\cdot(\mathbf{u}c) = \mathbf{u}\cdot\nabla c$. The discrete scheme below keeps the flux form so that the total amount of scalar is conserved exactly. This also enables us to design a correctness check for parallel computing in terms of volume preservation through time.
- There is no source term, to keep the problem simple.

With $A = 0$ the flow vanishes ($\mathbf{u} = 0$) and, with Dirichlet walls, the model reduces exactly to the heat equation with $\kappa$ as the thermal diffusivity. This special case provides an analytical reference.

All implementations solve the same discrete problem above.

## Grid and indexing

We use a cell-centred finite-volume grid. Each grid value is the average of $c$ over one cell, and the physical walls lie on cell faces.

- $n_x$: global cells along $x$ (columns); $n_y$: global cells along $y$ (rows).
- Step size $\Delta x = l_x/n_x$, $\Delta y = l_y/n_y$.
- Grid access uses $(i, j)$ = (row, column); top-down row ordering is used in field files.
- For global row $i=1,\dots,n_y$, the cell centre is $y = l_y-(i-\tfrac12)\Delta y$; for column $j=1,\dots,n_x$, it is $x = (j-\tfrac12)\Delta x$.
- Row $i=0$, row $i=n_y+1$, column $j=0$ and column $j=n_x+1$ are ghost cells outside the walls. The top wall is the face between rows 0 and 1, at $y = l_y$.
- The east face of cell $(i,j)$ is at $x = j\,\Delta x$. The north face of cell $(i,j)$ is at $y = l_y-(i-1)\Delta y$.

The five-point update never reads a corner ghost cell, so corner values need no convention. This replaces the earlier node-based grid with $\Delta x = l_x/(n_x+1)$ and its corner rule.

## Velocity on faces

Store or evaluate one normal velocity $(u,v)=(\frac{\partial\psi}{\partial y},\frac{\partial\psi}{\partial x})$ per cell face, computed from $\psi$ at the face's two end points (cell corners):

on east face, we compute,
$$
U_e(i,j) = \frac{\psi\big(j\Delta x,\; l_y-(i-1)\Delta y\big) - \psi\big(j\Delta x,\; l_y-i\Delta y\big)}{\Delta y}
$$

on thenorth face, we compute,
$$
V_n(i,j) = -\,\frac{\psi\big(j\Delta x,\; l_y-(i-1)\Delta y\big) - \psi\big((j-1)\Delta x,\; l_y-(i-1)\Delta y\big)}{\Delta x}
$$

The west and south faces of a cell are the east face of $(i,j-1)$ and the north face of $(i+1,j)$. With this construction the discrete divergence of every cell is exactly zero. We prove this in the following.
$$
(\nabla\cdot\mathbf u)_{i,j}
\approx
\frac{U_e-U_w}{\Delta x}
+
\frac{V_n-V_s}{\Delta y}
=\frac{1}{\Delta x \Delta y}[\psi_{NE}-\psi_{SE}
$$
$$
-\psi_{NW}+\psi_{SW}-\psi_{NE}+\psi_{NW}+\psi_{SE}-\psi_{SW}]=0
$$
For a closed container the streamfunction is zero on all four walls.

### Speed bounds

For the double gyre the velocity components satisfy, for all $t$,

$$
\lvert u\rvert \le u_{\max} = \frac{\lvert A\rvert\pi}{l_y},
\qquad
\lvert v\rvert \le v_{\max} = \frac{2\lvert A\rvert\pi\,(1+2\varepsilon)}{l_x}.
$$

Each face velocity is a difference quotient of $\psi$, so it is bounded by the maximum of the corresponding continuous velocity component. 

We can derive the bounds in the following. While $X=2x/l_x, \, Y=y/l_y$, for horizontal velocity, using $u=\partial \psi/ \partial y$, we have

$$u=\frac{A\pi}{l_y}\sin(\pi f)\cos(\pi Y)\leq \frac{|A|\pi}{l_y}
$$

For vertical velocity, using $v=-\partial \psi/\partial x$, we have 

$$v=-\frac{2A\pi}{l_x}\cos(\pi f)\sin(\pi Y)\frac{\partial f}{\partial X}\leq \frac{2|A|\pi(1+2\varepsilon)}{l_x}
$$

## Discretization

Write the update in conservative form. For each face, compute one total flux, advective plus diffusive $F_{\text{total}}=F_{\text{advection}}+F_{\text{diffusion}}$, using first-order upwinding for advection and a central difference for diffusion:

$$
F_e(i,j) = U_e \cdot \begin{cases} c_{i,j} & U_e > 0 \\ c_{i,j+1} & \text{otherwise} \end{cases}
\;-\; \kappa\,\frac{c_{i,j+1} - c_{i,j}}{\Delta x}
$$

$$
G_n(i,j) = V_n \cdot \begin{cases} c_{i,j} & V_n > 0 \\ c_{i-1,j} & \text{otherwise} \end{cases}
\;-\; \kappa\,\frac{c_{i-1,j} - c_{i,j}}{\Delta y}
$$

$$
\begin{align*}
c^{\,n+1}_{i,j} &= c^{\,n}_{i,j}
- \frac{\Delta t}{\Delta x}\big(F_e(i,j) 
- F_e(i,j-1)\big)\\
&- \frac{\Delta t}{\Delta y}\big(G_n(i,j) - G_n(i+1,j)\big)
\end{align*}
$$


Use forward Euler in time. The flow is evaluated at the old time level $t_n$ for time-dependent flows.

$$c^{n+1}=c^n+\Delta t\,\, \mathcal{L}(c^n, t_n)$$

Every face flux must be computed by the same expression, with the same operand order, in both cells that share the face, including cells on different ranks. 

Upwinding is stable and non-oscillatory but adds numerical diffusion of roughly $\lvert u\rvert\Delta x/2$ in $x$ and $\lvert v\rvert\Delta y/2$ in $y$. The leading spatial error is therefore first order in $\Delta x$ and $\Delta y$. Here is a detailed breakdown of the overall introduced error from the advection computation as well as other possible sources.

As we use first-order upwinding difference for advection term, we have 

$$\frac{\partial c}{\partial x}\approx \frac{c_{i,j}-c_{i, j-1}}{\Delta x}$$

where using a Taylor expansion, 
$$
\begin{align*}
\frac{c(x)-c(x-\Delta x)}{\Delta x}
&=\frac{c(x)-(c(x)-\Delta x c_x + c_{xx}\Delta x^2/2 +O(\Delta x^3))}{\Delta x}\\
&=c_x-\frac{\Delta x}{2}c_{xx}+O(\Delta x^2)
\end{align*}
$$
The leading spatial discretization error term here is $-\frac{\Delta x}{2}c_{xx}\approx O(\Delta x)$. For our two-dimensional case, this would introduce $O(\Delta x)+O(\Delta y)$.

Now, consider the diffusion, as we employ another discretization scheme, the central difference, we have 

$$
c_{xx}\approx
\frac{c_{i,j+1}-2c_{i,j}+c_{i,j-1}}{\Delta x^2}
$$

giving an error of $O(\Delta x^2)$, and a total error term of $O(\Delta x^2)+O(\Delta y^2)$ in two dimensions. 

Finally, we have the error introduced from the time integration, giving $O(\Delta t)$. Indeed, we can reduce this amount of error by introducing a higher-order discretization scheme, however, that would require higher computational capacity and is thus left to our future study if possible. 

## Boundary conditions

Each run adopts the following rule for all four walls.

`noflux`: no scalar crosses a wall. Before each step, set every ghost cell to the value of its adjacent interior cell, $c_{\text{ghost}} = c_{\text{adjacent}}$, so $F_{wall}^{diff}=-\kappa \frac{c_i-c_g}{\Delta x}=0$. With zero wall velocity, we have $F_{wall}^{adv}=0$. Togher, the total wall flux is zero.

As ghost cells on physical walls now depend on the current interior values. They must be **refilled at the start of every step**, not only at initialization. 

Note that we don't actually need to take care of the four corner ghost cells as they are not involved in any step of the computation.

## Initial conditions

Initialize interiors from cell-centre global coordinates so that a change in decomposition cannot change the initial condition.

| `case` | Initial field | Notes |
| --- | --- | --- |
| `heated_top` | $c_0$ everywhere | Heat-only demonstration. flow amplitude $A = 0$, Dirichlet walls with a hot top. |
| `constant` | $c_0$ everywhere | With all no-flux walls equal $c_0$. Must stay constant under the flow, even when the double gyre is active. |
| `sine_decay` | $A_0\sin\!\big(\tfrac{\pi x}{l_x}\big)\sin\!\big(\tfrac{\pi y}{l_y}\big)$ | Flow amplitude $A = 0$. Analytical reference. |
| `pigment_drop` | $c_0 + A_0\exp\!\Big(-\dfrac{(x-x_0)^2 + (y-y_0)^2}{2\sigma^2}\Big)$ | Default advection–diffusion case. A Gaussian-shaped pigment blob centered at $(x_0,y_0)$, with $\sigma$ controlling its width. The double gyre will transport, stretch this blob throughout time. |

## Time step
This part, we explain how to choose an appropriate time step. As we employ the forward Euler scheme, we must constrain the time step size for stability.
Define the advective and diffusive rates:

$$
a = \frac{u_{\max}}{\Delta x} + \frac{v_{\max}}{\Delta y}
\quad\text{($u_{\max}$, $v_{\max}$ from the flow's speed bounds)},
\qquad
d = 2\kappa\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right)
$$

Validate positive dimensions, positive $\kappa$, positive $\Delta t$ and

$$
\Delta t\,(a + d) \le 1
$$

This serves as a sufficient stability restriction given our schemes. For a discretely divergence-free flow this keeps every update weight non-negative in cells away from Dirichlet walls. Each new value is then a convex combination of old values, so no new maximum or minimum is created there and the scheme is stable. With $\mathbf{u} = 0$ the condition reduces to the earlier $r_x + r_y \le \tfrac12$, where $r_x = \kappa\Delta t/\Delta x^2$ and $r_y = \kappa\Delta t/\Delta y^2$.

For automatic time stepping, `cfl` in $(0,1]$ is a safety fraction:

$$
\Delta t = \frac{\mathrm{cfl}}{a + d}
$$

We can then try different time step size to investigate the stability dynamics.

Record the Péclet numbers for each run, $\mathrm{Pe} = UL/\kappa$, using a characteristic speed $U$ and length $L$, and the grid Péclet number $\mathrm{Pe}_{\text{grid}} = \max(u_{\max}\Delta x,\, v_{\max}\Delta y)/\kappa$. A large Péclet number means the transport by the advection dominates diffusion. When $\mathrm{Pe}_{\text{grid}}$ is well above 2, upwind numerical diffusion is comparable to or larger than the physical diffusion, and this means while we can keep the computation stable, it may become physically inaccurate already.

This constraint is derived in the following way. As we use first-order upwind difference scheme for advection, we have numerical diffusion of roughly $\kappa_{num}\approx \frac{|u|\Delta x}{2}$, and similarly for $y$. Then, we have $\kappa_{num}/\kappa \approx \frac{|u|\Delta x}{2\kappa}$.

Finally, we can represent the pipeline for fitting a proper time step size. Given parameters $\Delta x, \Delta y, \kappa, A, \varepsilon, l_x, l_y$, we can compute speed bounds $u_{max}, v_{max}$ and the following rates $a, d$ defined above. As a last step, select $\Delta t=\text{cfl}/(a+d)$ and put it into test.

## Data storage and shared computation

Use contiguous `std::vector<double>` storage. Each local field includes a one-cell halo on every side, and each solver owns a current and next field. The row stride is `local_nx+2`.

For steady flows, store local face velocities: `U` has `local_ny` rows of `local_nx+1` east/west faces and `V` has `local_ny+1` rows of `local_nx` north/south faces, all computed once from global coordinates. For time-dependent flows, recompute the local face velocities at the start of each step; this cost is part of the measured solver loop. Record which approach is used, because it changes the memory traffic and arithmetic per cell.

The shared stencil computes a single destination value from the five cell values and the four face velocities of that cell. It contains no MPI calls, OpenMP regions, allocation or output. The serial and parallel solvers control the loops. This separation makes formula changes explicit and keeps the serial executable free from MPI/OpenMP dependencies.

Initialize interiors by global coordinates. Fill physical-wall ghosts by the boundary rule and process halos by exchange before their values are used in each step.

## Domain decomposition

### Rows

Distribute $n_y$ rows over the $p$ ranks. With $q = \lfloor n_y/p\rfloor$ and $r = n_y \bmod p$, the first $r$ logical row positions get $q+1$ rows and the rest get $q$. Compute global starting rows from prefix sums or the equivalent arithmetic formula.

Use a process topology with upper/lower neighbors. Do not confuse Cartesian rank numbers with world rank numbers when reordering is enabled.

### Two-dimensional blocks

Use a $[p_y, p_x]$ Cartesian process grid, with $p_y p_x = p$. Divide rows and columns independently, including remainders. An explicit grid must fit within the global grid dimensions; reject zero-sized local blocks. Automatic process-grid selection should also respect this requirement.

Derive neighbors with MPI topology utilities. Physical walls belong only to ranks on the global edge. Columns are not contiguous in row-major storage: implement either column packing or an appropriate MPI datatype, document it and keep it consistent across communication-mode comparisons.

## Halo exchange

Upwinding and the five-point diffusion term use only the four nearest neighbours, so a one-cell halo is sufficient.

The blocking baseline can use `MPI_Sendrecv` with consistent partners and tags. It must not rely on eager buffering to avoid deadlock. Sides on a physical wall have no rank neighbour; fill those ghosts by the boundary rule instead.

For overlap:

1. Post receives into the current field's halos.
2. Post sends from its owned edge values or dedicated send buffers.
3. Fill physical-wall ghosts from the current owned values; this is local work.
4. Update points that do not depend on incoming halos.
5. Complete all sends and receives.
6. Update the remaining owned points exactly once.
7. Swap current and next buffers.

Do not modify send buffers, read incomplete receive buffers or reuse communication resources before the matching requests complete. Small local blocks may have no strict interior; handle them without duplicate edge updates.

Nonblocking calls do not guarantee a performance improvement. MPI progress, message size and available independent work determine whether overlap helps.

## OpenMP and MPI thread support

Begin with `MPI_Init_thread(..., MPI_THREAD_FUNNELED, ...)` and verify support. Only the initializing thread performs MPI calls. OpenMP worksharing distributes destination grid points among threads. Each thread computes the face fluxes it needs for its own destination cells; do not write shared face-flux arrays from several threads. Complete the worksharing region before the master thread swaps buffers or starts dependent work.

Support builds without OpenMP by guarding OpenMP-specific calls with the build definition. Do not add global per-step barriers without a demonstrated need; halo completion and local worksharing must already enforce their own dependencies.

## Diagnostics

Total scalar over interior cells:

$$
M = \Delta x\,\Delta y \sum_{i,j} c_{i,j}.
$$

With `noflux` walls $M$ is conserved up to rounding. In the parallel solver $M$ needs an `MPI_Allreduce`; compute it only at snapshot intervals or after the timed loop, and compare serial and parallel values with a tolerance because the summation order differs.

## Output and metrics

Small validation runs may gather the global matrix. Gathering and file writing happen outside solver-loop timing. Global output contains the $n_y \times n_x$ interior cells once, never ghost cells or duplicate process halos. Record the wall type and wall values in the run metadata.

Time the whole loop, including necessary communication, wall-ghost filling and any per-step velocity evaluation. Synchronize before the timed MPI region and report the maximum elapsed time over ranks. Keep initialization, validation and output timings separate. Phase timings under overlap are diagnostic and must not be blindly added together.

## Scope

The project uses CPUs, a five-point stencil with first-order upwind advection, constant diffusivity, a prescribed divergence-free velocity field and rectangular geometry. Solving for the velocity (Navier–Stokes or Boussinesq flow), higher-order advection schemes with wider halos, variable materials, complex geometry, implicit methods, CUDA and a web interface are outside the current plan.

## Follow-up from the advection–diffusion change

This design supersedes the heat-only, node-based design. The following files still describe the old model and need matching updates:

- `docs/experiment_plan.md`: new checks for mass conservation and constant field under the double-gyre flow; the sine-decay test ($A = 0$) now uses cell centres.
- `config/defaults.json`, `config/README.md`: `kappa`, the double-gyre parameters ($A$, $\varepsilon$, $\omega$), `boundary_type`, `pigment_drop` parameters, the new time-step rule and $\Delta x = l_x/n_x$.
- `experiments/*.json`: flow and wall settings; weak scaling now uses $l_x = p_x\,n_{x,\text{local}}\,\Delta x$, and needs fixed speed bounds as the domain grows.
- `src/*/README.md`, `results/raw/metrics_columns.csv`, `scripts/README.md`, `hpc/job_template.sh`, `report/report.md` and the top-level `README.md`.
