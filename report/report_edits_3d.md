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

Figure numbers after the change: Figure 1 the 2D domain (redrawn as one panel), Figure 2 the 3D domain, Figure 3 the 2D flow evolution, Figure 4 the 3D evolution.

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

In 3D the streamfunction and vorticity are no longer scalars, so we solve for the velocity
$\mathbf{u}=(u,v,w)$ and the pressure $\phi$ (divided by the constant density) in the basin
$\Omega=[0,L_x]\times[0,L_y]\times[0,L_z]$:
\begin{equation}
\frac{\partial \mathbf{u}}{\partial t}+\nabla\cdot(\mathbf{u}\mathbf{u})
=-\nabla\phi+\nu\nabla^2\mathbf{u}+\mathbf{f},\qquad
\nabla\cdot\mathbf{u}=0,\qquad
\frac{\partial c}{\partial t}+\nabla\cdot(\mathbf{u}c)=\kappa\nabla^2 c .
\label{eq:ns3d}
\end{equation}
The walls are impermeable and free-slip, the pigment is passive with $\partial c/\partial n=0$,
and the wind pushes in $x$:
\begin{equation}
f_x=\frac{F_0}{k_y}\cos(k_y y)\,g(z),\qquad k_y=\frac{2\pi}{L_y},\qquad
g(z)=\frac{4\,e^{-4(1-z/L_z)}}{1-e^{-4}} .
\label{eq:wind3d}
\end{equation}
Its curl is $F_0\sin(k_y y)\,g(z)$, the 2D forcing weighted by a profile of mean one that
concentrates it near the surface (Figure~\ref{fig:domain-3d}). The flow therefore varies with
depth, and the pressure creates a vertical velocity.

The pigment and $\phi$ live at the cell centres and $u$, $v$, $w$ at the faces (a staggered
grid). Each step updates the pigment with the scheme of Section~3.1, now with a seven-point
stencil, predicts a velocity $\mathbf{u}^*$ with the Adams--Bashforth scheme, and projects it:
\begin{equation}
A\phi=-\frac{D\mathbf{u}^*}{\Delta t},\qquad A=-DG,\qquad
\frac{\partial\phi}{\partial n}=0,\qquad
\mathbf{u}^{n+1}=\mathbf{u}^*-\Delta t\,G\phi,
\label{eq:projection}
\end{equation}
where $D$ and $G$ are the discrete divergence and gradient, so that $D\mathbf{u}^{n+1}=0$.
This is the 3D Poisson problem of the project, solved by the CG of Section~3.2 with the
seven-point Laplacian. Under the Neumann condition $A$ is singular, as constants lie in its
null space, so we remove the mean of the right-hand side, which is zero in exact arithmetic,
and keep $\phi$ at zero mean, where CG applies.

In parallel, the grid is split into blocks in all three directions, with a one-cell halo on each
face. A CG iteration thus exchanges six faces instead of four edges, and still needs two global
sums. Exchanges are blocking, and threads share the loops within a rank as in 2D. The model is
a separate program, \texttt{model\_b}, with no GPU version.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.45\linewidth]{figures/domain_b.pdf}
\caption{Domain and initial condition of the three-dimensional benchmark: the basin
$\Omega=[0,1]^3$ (surface $z=1$ shaded) and the initial pigment $c(x,y,z,0)$, drawn as points
whose density follows $c$, a Gaussian blob of amplitude $A_0=1$ and width $\sigma=0.08$ centred
at $(x_0,y_0,z_0)=(0.3,0.5,0.5)$. The dotted curves are the two factors of the wind force
$f_x=(F_0/k_y)\cos(k_y y)\,g(z)$: $\cos(k_y y)$ on the floor, about $x=0.5$, and the depth profile
$g(z)$ on the wall $x=0$. All walls are impermeable and free-slip, with no pigment flux,
$\partial c/\partial n=0$ for the outward normal $\mathbf{n}$; the flow starts from rest.}
\label{fig:domain-3d}
\end{figure}
```

---

## 4.2 Flow and Pigment Evolution (replace the whole subsection)

Both evolution figures sit here: your Figure 3 (2D, unchanged) and the 3D one (`domain_3d.pdf`, now Figure 4).

```latex
\subsection{Flow and Pigment Evolution}

The checks of Section~4.1 compare numbers on small grids. Here we check that one full run
behaves as the equations say it should, in 2D and in 3D. Both runs use $F_0=0.8$,
$\nu=0.004$ and $\kappa=10^{-4}$, and a spin-up of 8 time units before the blob is released;
$\tau$ is the time since the release.

In 2D (Figure~\ref{fig:flow}, $128\times128$ cells), the wind changes sign at mid-height and
drives two gyres of opposite sense, which both carry fluid to the left between them, where the
blob sits. The pigment reaches the left wall, splits into one branch per gyre, and by
$\tau=9.77$ has travelled around both. Its peak falls from 1 to 0.06 while its mass is constant
to rounding: advection along the streamlines, diffusion across them, and no loss through the
walls.

In 3D (Figure~\ref{fig:flow-3d}, $48\times48\times24$ cells, depth 0.5), the pigment takes the
same path seen from above, and the depth adds two effects. First, the wind acts near the
surface, so the flow is faster there (0.53 against 0.14 near the floor, at the release point):
by $\tau=2.5$ the upper pigment has already turned back along the walls, while the lower
pigment is still at the left wall. Second, the surface current sinks where it meets the left
wall, with a vertical velocity of up to 0.37, and carries 74\% of the pigment into the lower
half. The mass is again conserved, to a relative $2\times10^{-14}$, and the divergence stays
below $5\times10^{-12}$ after every projection.

% Figure 3 (2D flow evolution): unchanged

\begin{figure}[htbp]
\centering
\includegraphics[width=\linewidth]{figures/domain_3d.pdf}
\caption{Flow and pigment in the three-dimensional model on a $48\times48\times24$ grid in
$[0,1]\times[0,1]\times[0,0.5]$, with $F_0=0.8$, $\nu=0.004$, $\kappa=10^{-4}$ and a spin-up of
8 time units. (a) Streamlines of the depth-averaged horizontal flow at the release, on the
surface, with the blob of Figure~\ref{fig:domain-3d}. (b)--(d) Pigment concentration $c$ at times
$\tau$ after the release, as points whose density follows $c$, on a common square-root colour
scale.}
\label{fig:flow-3d}
\end{figure}
```

---

## Figure 1 (replace the figure; the panel (b) is now inside the basin)

```latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.45\linewidth]{figures/domain.pdf}
\caption{Domain and initial condition of the fixed-size benchmark: the basin
$\Omega=[0,1]^2$ and the pigment $c(x,y,0)$, a Gaussian blob of amplitude $A_0$ and width
$\sigma$ centred at $(x_0,y_0)$, drawn as points whose density follows $c$. The dotted curve is
the wind profile $\sin(2\pi y/L_y)$ of the forcing $F=F_0\sin(2\pi y/L_y)$, which changes sign at
$y=L_y/2$. On the walls $\partial\Omega$, $\psi=\zeta=0$ (impermeable, free-slip) and
$\partial c/\partial n=0$.}
\label{fig:domain}
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
