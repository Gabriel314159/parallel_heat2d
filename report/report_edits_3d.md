# Report edits for the three-dimensional model (Model B)

These blocks replace everything about the old 3D extension (Model A, the pigment in the 2D flow) with the real 3D model, which has its own velocity and a pressure projection. Paste each block over the matching part of the `.tex` file. Labels used here:

| Label | What it points to |
| --- | --- |
| `sec:model3d` | 3.4 |
| `eq:ns3d`, `eq:wind3d`, `eq:projection` | the equations of 3.4 |
| `fig:domain-3d` | new Figure 2, the 3D overview (`results/figures/domain_b.pdf`) |
| `fig:flow-3d` | the 3D evolution figure, now in 4.2 (`results/figures/domain_3d.pdf`) |
| `sec:setup-3d`, `tab:setup-3d` | 5.3.8 and its table (your Table 5) |
| `tab:strong-3d`, `fig:strong-3d` | the 3D strong-scaling table and figure in 6.9 |
| `fig:domain`, `fig:flow` | your Figures 1 and 3 (the 2D ones); match these to your labels |

Figure numbers after the change: Figure 1 the 2D domain, Figure 2 the 3D domain, Figure 3 the 2D flow evolution, Figure 4 the 3D evolution.

Files made for these edits:

- `results/figures/domain_b.pdf`, from `scripts/plot_domain_b.py`: the new Figure 2.
- `results/figures/domain_3d.pdf`, from `scripts/plot_domain_3d.py`, now redrawn from Model B runs: the 3D evolution figure for 4.2.
- `results/figures/strong_scaling_3d_table.tex`, from `scripts/make_table_strong_3d.py`: the table for 6.9.
- `results/figures/strong_scaling_3d.pdf` (already there): the figure for 6.9.

---

## Section 3 opening (one sentence)

In the paragraph that lists the subsections, replace "and Section 3.4 the three-dimensional extension of the pigment" with:

```latex
and Section~\ref{sec:model3d} the three-dimensional model, whose flow is computed with a pressure projection.
```

---

## 3.4 Three-Dimensional Model with Pressure Projection (replace the whole subsection)

```latex
\subsection{Three-Dimensional Model with Pressure Projection}
\label{sec:model3d}

In three dimensions the streamfunction and the vorticity are no longer scalars, so the
formulation of Section~2 does not carry over. We therefore solve the 3D model in the
velocity $\mathbf{u}=(u,v,w)$ and a pressure $\phi$ (pressure divided by the constant
density). In the basin $\Omega=[0,L_x]\times[0,L_y]\times[0,L_z]$,
\begin{equation}
\frac{\partial \mathbf{u}}{\partial t}+\nabla\cdot(\mathbf{u}\mathbf{u})
=-\nabla\phi+\nu\nabla^2\mathbf{u}+\mathbf{f},\qquad
\nabla\cdot\mathbf{u}=0,\qquad
\frac{\partial c}{\partial t}+\nabla\cdot(\mathbf{u}c)=\kappa\nabla^2 c .
\label{eq:ns3d}
\end{equation}
The pigment is passive, as in 2D. All six walls are impermeable and free-slip: the normal
velocity is zero and the tangential velocity has zero normal derivative; the pigment satisfies
$\partial c/\partial n=0$, so its total mass is conserved. The wind pushes in $x$,
\begin{equation}
f_x=\frac{F_0}{k_y}\cos(k_y y)\,g(z),\qquad k_y=\frac{2\pi}{L_y},\qquad
g(z)=\frac{4\,e^{-4(1-z/L_z)}}{1-e^{-4}},
\label{eq:wind3d}
\end{equation}
with $f_y=f_z=0$. Its curl, $-\partial f_x/\partial y=F_0\sin(k_y y)\,g(z)$, is the forcing of
Section~2 weighted by the profile $g$, whose mean over the depth is one: the 2D wind,
concentrated near the surface. Because the forcing changes with depth, the horizontal flow does too,
and the pressure creates a vertical velocity that no force imposes directly.
Figure~\ref{fig:domain-3d} summarizes the domain, the initial pigment and the forcing.

The grid has $n_x\times n_y\times n_z$ cells. The pigment and $\phi$ live at the cell centres,
and $u$, $v$, $w$ at the centres of the faces normal to $x$, $y$, $z$ (a staggered, or MAC,
grid). One time step has four parts.
\begin{enumerate}
\item The pigment is updated with the current velocity, by the scheme of Section~3.1 with
fluxes through the six faces of a cell, a seven-point stencil.
\item A predicted velocity $\mathbf{u}^*$ is computed from the advection, viscosity and wind,
with the Adams--Bashforth scheme of Section~3.1. It does not yet satisfy $\nabla\cdot\mathbf{u}=0$.
\item The pressure solves
\begin{equation}
A\phi=b,\qquad A=-DG,\qquad b=-\frac{D\mathbf{u}^*}{\Delta t},\qquad
\frac{\partial\phi}{\partial n}=0 \ \text{on}\ \partial\Omega,
\label{eq:projection}
\end{equation}
where $D$ and $G$ are the discrete divergence and gradient of the staggered grid, so that
$-DG$ is the seven-point Laplacian.
\item The velocity is corrected, $\mathbf{u}^{n+1}=\mathbf{u}^*-\Delta t\,G\phi$, which makes
$D\mathbf{u}^{n+1}=0$ up to the tolerance of the solve.
\end{enumerate}

Equation~\eqref{eq:projection} is the three-dimensional Poisson problem of this project. It is
solved by the matrix-free CG of Section~3.2, with the seven-point stencil in place of the
five-point one and the previous $\phi$ as the initial guess. One difference matters: with the
Neumann condition, $A$ is singular, since adding a constant to $\phi$ does not change $A\phi$.
A solution exists only if the mean of $b$ is zero, which holds because no fluid crosses the
walls. We remove the mean of $b$, fix the mean of $\phi$ to zero, and run CG in the space of
mean-zero vectors, where $A$ is positive definite. The stopping criterion is that of
Equation~(20), with $b$ in place of $\zeta^{n+1}$.

In parallel, the grid is split in all three directions into blocks (slabs and pencils are also
available), and the process grid is chosen to minimize the area of the faces between blocks.
Each block keeps one layer of halo cells on each of its six faces. In every CG iteration the
search direction exchanges these six faces, instead of the four edges of 2D, and the two dot
products need two global sums, as in 2D. The velocity halos are exchanged one direction after
another, so that the cells along the edges of a block, which the momentum advection needs,
arrive as well, as the corner cells do in 2D. Exchanges are blocking, and within a rank OpenMP
threads share the loops. The model is a separate program, \texttt{model\_b}, with a serial build
and an MPI + OpenMP build; it has no GPU version.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.6\linewidth]{figures/domain_b.pdf}
\caption{Computational domain and initial conditions of the three-dimensional benchmark.
(a) The basin $\Omega=[0,1]^3$, with the surface $z=1$ shaded, and the initial pigment
concentration $c(x,y,z,0)$, drawn as a cloud of points whose density follows $c$ and coloured
on the scale of Figure~\ref{fig:domain}. The pigment is a Gaussian blob of amplitude $A_0=1$
and width $\sigma=0.08$ centred at $(x_0,y_0,z_0)=(0.3,\,0.5,\,0.5)$. On the walls
$\partial\Omega$ the normal velocity $\mathbf{u}\cdot\mathbf{n}$ is zero, the tangential velocity
has zero normal derivative (free slip) and the pigment flux $\partial c/\partial n$ is zero,
where $\mathbf{n}$ is the outward normal. The flow starts from rest. The dotted curves show the two factors of the wind force
$f_x=(F_0/k_y)\cos(k_y y)\,g(z)$ of Equation~\eqref{eq:wind3d}, each measured from a thin
reference line: on the floor, $\cos(k_y y)$ against $y$, which changes sign with $y$ as the
forcing of Figure~\ref{fig:domain} does, pushing towards $-x$ at mid-height and towards $+x$
near the walls $y=0$ and $y=L_y$; on the back wall $x=0$, the depth profile $g(z)$, which
concentrates the wind near the surface and whose mean over the depth is one.}
\label{fig:domain-3d}
\end{figure}
```

---

## 4.2 Flow and Pigment Evolution (add after the last paragraph, with the figure)

The current Figure 2 (the 3D clouds) moves here. Its run uses the same wind, viscosity, spin-up and diffusivity as the 2D run, so the two figures can be read side by side.

```latex
The same check applies to the three-dimensional model. Figure~\ref{fig:flow-3d} shows a run on
a $48\times48\times24$ grid in a basin of depth $L_z=0.5$, with the parameters of
Figure~\ref{fig:flow} and the blob of Figure~\ref{fig:domain-3d}, centred at mid-depth. Seen from
above, the pigment follows the path of the 2D run: the flow between the gyres carries it to the
left wall, where it splits into two branches that turn along the top and bottom walls. The depth
adds what the 2D model cannot show. Since the wind acts mostly near the surface, the flow there
is faster: at the release point it moves left at 0.53 near the surface and at 0.14 near the
floor. By $\tau=2.5$ the pigment near the surface has therefore already travelled along the
walls and back towards the middle of the basin (its mean $x$ is 0.43 in the top layer), while
near the floor it is still at the left wall (mean $x$ of 0.10). The pigment also sinks: where
the surface current meets the left wall, the fluid has nowhere to go but down, and the vertical
velocity this creates, up to 0.37 and imposed by no force, carries 74\% of the pigment into the
lower half of the basin. By $\tau=10$ it has spread through the whole basin. Meanwhile the peak
concentration falls from 0.96 at the release to 0.14 and then to 0.026, the total pigment stays
the same to within a relative $2\times10^{-14}$, and the divergence of the velocity stays below
$5\times10^{-12}$ after every projection.

\begin{figure}[htbp]
\centering
\includegraphics[width=\linewidth]{figures/domain_3d.pdf}
\caption{Wind-driven flow and pigment transport in the three-dimensional model, on a
$48\times48\times24$ grid in the basin $[0,1]\times[0,1]\times[0,0.5]$, with $F_0=0.8$,
$\nu=0.004$, $\kappa=10^{-4}$ and a spin-up of 8 time units. The pigment concentration $c$ is
drawn as a cloud of points whose density follows $c$, on the colour scale of
Figure~\ref{fig:flow}, at times $\tau$ after the release; the shaded top face is the surface.
(a) The blob of Figure~\ref{fig:domain-3d}, centred at mid-depth $z=0.25$, at the release.
(b) The pigment carried to the left wall and along the side walls, faster near the surface, and
sinking at the left wall. (c) The pigment spread through the whole basin.}
\label{fig:flow-3d}
\end{figure}
```

---

## 5.3.8 Three-Dimensional Experiments (replace the whole subsection, including Table 5)

```latex
\subsubsection{Three-Dimensional Experiments}
\label{sec:setup-3d}

The three-dimensional experiments use \texttt{model\_b} (Section~\ref{sec:model3d}). Every job
first runs 69 correctness checks on its own build and times nothing unless all of them pass.
They compare the parallel fields $c$, $\phi$, $u$, $v$, $w$ with the serial ones in eight
layouts, and test the operators against known answers: a Neumann problem with a manufactured
solution, the spatial order of the pressure and momentum operators, the divergence left by
the projection, a constant pigment, a decaying cosine mode, and a uniform wind, under which the
flow must not depend on depth.

The fixed problem $P_3$ is the unit cube with $1024^3=1.07\times10^9$ cells and one time step
from rest, with the surface wind, $F_0=0.8$, $\nu=10^{-3}$ and $\kappa=10^{-4}$. Its single
pressure solve starts from zero and takes 3535 CG iterations. We use a tolerance of $10^{-8}$,
because $10^{-10}$ is at the limit of rounding for a system of this size. Strong scaling runs
$P_3$ on 1 to 64 ranks on one node and on 128 and 256 ranks on two and four nodes, one thread
per rank, with the placement of Section~\ref{sec:setup-strong-mpi}: up to 32 ranks on one
socket, then 32 per socket. Weak scaling keeps about $162^3\approx4.3\times10^6$ cells per
rank, from $162^3$ on one rank to $1024^3$ on 256, for two time steps. Table~\ref{tab:setup-3d}
lists the configurations.

\begin{table}[htbp]
\centering
\small
\begin{tabular}{p{0.2\linewidth}p{0.72\linewidth}}
\toprule
Item & Configuration \\
\midrule
Model & 3D velocity and pressure with projection, passive pigment (Section~\ref{sec:model3d}) \\
Validation & 69 checks per job; serial against 1 to 6 ranks with 1 or 2 threads, in slabs,
pencils and blocks \\
Fixed problem $P_3$ & $1024^3$ cells on the unit cube, one step from rest, CG tolerance $10^{-8}$ \\
Strong scaling & $P_3$ on 1, 2, 4, \dots, 256 ranks with one thread each, blocks, blocking
exchange; three repetitions at 64, 128 and 256 ranks, one run otherwise \\
Weak scaling & $n^3$ cells on $p$ ranks with $n=162$, 204, 256, 322, 406, 512, 646, 812, 1024
for $p=1$, 2, 4, \dots, 256; two steps, one warm-up and two repetitions \\
Status & 1 to 64 ranks done; 128 and 256 ranks and weak scaling queued on Vanda \\
\bottomrule
\end{tabular}
\caption{Configurations of the three-dimensional experiments.}
\label{tab:setup-3d}
\end{table}
```

---

## 6.9 Three-Dimensional Model and the 2D/3D Poisson Comparison (replace the whole subsection)

The table is `results/figures/strong_scaling_3d_table.tex`; paste it, or `\input` it, where marked. The 128 and 256 rank rows and the weak-scaling part will be added when job 1567813 finishes.

```latex
\subsection{Three-Dimensional Model and the 2D/3D Poisson Comparison}

Finally, we come to the three-dimensional model of Section~\ref{sec:model3d}. As with the 2D
solver, we first need to ensure that the parallel program solves the same problem as the serial
one. Every job passed all 69 checks of Section~\ref{sec:setup-3d} before it timed anything.
The parallel fields agree with the serial ones to within $1.2\times10^{-13}$ in all eight
layouts, the pressure and momentum operators converge at the expected second order (measured
orders between 2.01 and 2.05), and after each projection the divergence of the velocity is at
most $1.4\times10^{-12}$.

% results/figures/strong_scaling_3d_table.tex

\begin{figure}[htbp]
\centering
\includegraphics[width=\linewidth]{figures/strong_scaling_3d.pdf}
\caption{Strong scaling of the three-dimensional model over MPI ranks on $P_3$, on one node:
runtime, speedup and efficiency. The hollow point is the extrapolated one-rank run.}
\label{fig:strong-3d}
\end{figure}

Table~\ref{tab:strong-3d} shows the strong scaling of $P_3$ within one node. As in 2D, the
Poisson solve dominates the loop, with $T_P/T$ at least 99.3\%, and every run takes the same
3535 CG iterations, so $T$ and $t_{\mathrm{it}}$ tell the same story. The scale of the problem
is, however, of a different order: a single step takes more than eleven hours on one core
and still 20 minutes on 64 ranks.

Figure~\ref{fig:strong-3d} shows that the runtime follows the ideal curve closely up to 8 ranks,
where the efficiency is still 0.94. What differs from 2D is the role of communication. At 64
ranks the halo exchanges and global sums take 15.4 and 13.8 ms of an iteration of 330 ms, about
9\%, whereas in 2D they took 40\% (Section~6.3). This is because each rank now owns a block of
$256^3\approx1.7\times10^7$ cells, whose update takes about 300 ms per iteration, which is long
enough for the messages and the sums to be cheap in comparison.

The efficiency nevertheless drops to 0.80 at 16 ranks and 0.55 at 32, and the cause has to be
the computation itself. If we take the halo and sum out of $t_{\mathrm{it}}$, the remaining time
per iteration at 32 ranks is 616 ms, 1.8 times the ideal value scaled from 2 ranks. Recall that
up to 32 ranks all of them sit on one socket, while 64 ranks use both. Doubling the ranks from
16 to 32 within one socket gives a speedup of only 1.38, but doubling them from 32 to 64,
which also doubles the sockets, gives 2.0, and the efficiency stays at 0.55. This suggests that
the 3D solve is limited by the memory bandwidth of a socket: one field of $1024^3$ values takes
8.6 GB, so no block fits in cache at any of these rank counts, and, unlike the 2D benchmark, no
superlinear point appears. As in Section~6.4, we have not confirmed this with hardware
counters.

The runs on two and four nodes and the weak-scaling ladder of Table~\ref{tab:setup-3d} are
still in the queue on Vanda, and their results will complete this section.

\paragraph{2D against 3D Poisson solves.} The pressure equation~\eqref{eq:projection} also
gives us the comparison of a 2D and a 3D Poisson solver that Lecture~12 asks for, at the same
resolution of 1024 cells in each direction. Three differences stand out.
\begin{enumerate}
\item One CG iteration costs 10.8 ns per cell on one core in 3D, against 2.4 ns in 2D. The
seven-point stencil explains only a small part of this. The rest is most likely memory: a 2D
CG vector of $1024^2$ values (8 MB) fits in the cache of a socket, while a 3D one (8.6 GB) has
to stream from main memory, as on the $4096^2$ grid of Section~6.8.
\item The 3D solve needs 3535 iterations, against 1095 per step in 2D. The two numbers are not
directly comparable, since the 3D solve starts from zero with a tolerance of $10^{-8}$ and the
2D ones from the previous $\psi$ with $10^{-10}$. In both cases, however, the condition number
of the Laplacian grows like $n^2$ on a grid of $n$ cells per direction, so the iteration count
of unpreconditioned CG grows like $n$, and a larger grid costs more iterations in 3D and 2D
alike.
\item A 3D block exchanges six faces per iteration rather than four edges, about
$3.9\times10^5$ values for a $256^3$ block, while the global sums stay at two per iteration.
The communication still takes a smaller share of an iteration than in 2D, because the work of a
block grows with its volume and its halo only with its surface.
\end{enumerate}
```

---

## Other places that still describe the old 3D extension

Not changed here; a sentence each will do once the remaining results are in.

- **Opening of Section 6:** "before closing with the GPU and the three-dimensional extension" could say "the three-dimensional model".
- **Section 7, Numerical limitations:** "the 3D extension keeps a two-dimensional flow with no vertical velocity" is no longer true. A replacement: "The 3D model uses free-slip walls and no buoyancy, so the pigment does not act on the flow, and its pressure solve, like the 2D one, has no preconditioner."
- **Section 7, Scope:** "No 3D timings and no 3D Poisson solver are included" should go. Say instead that the 3D runs on more than one node and the 3D weak scaling are missing (until job 1567813 finishes).
- **Section 8, Conclusion:** add one sentence on the 3D strong scaling (speedup 35.1 at 64 ranks, efficiency 0.55, limited by memory bandwidth rather than communication).
- **Section 4.1** describes only the 2D checks. The 3D checks are now in 5.3.8 and their results in 6.9. If you prefer all checks in 4.1, move the first paragraph of 5.3.8 there.
