# Report edits for Sections 4, 5 and 6

Each block below is LaTeX to paste over the matching part of the working `.tex` file. I worked from the PDF, not from your source, so three things need adapting when you paste:

- **Labels.** I used the labels in the table below. Rename them to yours, or add them.
- **Tables and figures I did not change** are marked `% unchanged` with their caption title; keep your existing code there.
- **Notation.** Section 5 used $P$, $T$, $C$ for processes, threads and cores, while Section 6 uses $p$ for the core count and $T$ for the loop time. I rewrote Section 5 in the notation of Section 6, so the definitions that moved out of 6.2 now sit in 5.2.

| Label | What it points to |
| --- | --- |
| `sec:validation` | 4.1 Validation Strategy |
| `sec:environment`, `sec:measurement` | 5.1, 5.2 |
| `sec:setup-serial`, `sec:setup-strong-mpi`, `sec:setup-openmp`, `sec:setup-weak`, `sec:setup-hybrid`, `sec:setup-comm`, `sec:setup-gpu` | 5.3.1 to 5.3.7 |
| `sec:res-serial`, `sec:res-hybrid` | 6.1, 6.5 |
| `tab:environment`, `tab:benchmark` | Tables 2 and 3 (unchanged) |
| `tab:layouts`, `tab:comm-configs` | Tables 4 and 5 (unchanged) |
| `tab:validation` | Table 1, moved to 6.0 |
| `tab:strong-mpi`, `tab:strong-openmp`, `tab:weak-mpi`, `tab:hybrid-layouts`, `tab:hybrid-node`, `tab:node-spread` | as in `results/figures/*.tex` |
| `tab:placement`, `tab:comm`, `tab:gpu` | Tables 13, 17 and 18 (unchanged) |
| `fig:flow`, `fig:weak-breakdown` | Figure 2 (new file), the new weak-scaling figure |
| `eq:speedup`, `eq:weak-factor` | Equations (22) and (27) |

Files made for these edits:

- `results/figures/flow_evolution.pdf` and `.png`: the new Figure 2, from `scripts/plot_flow_evolution.py`.
- `results/figures/weak_scaling_mpi_table.tex`: Table 12 with nine columns instead of twelve, from `scripts/make_table_weak_mpi.py`.

---

## Section 4

### 4.1 Validation Strategy (replace the whole subsection)

The strategy now lives here and the numbers live in 6.0. Table 1 moves to 6.0.

```latex
\subsection{Validation Strategy}
\label{sec:validation}

Before any timing, we need to ensure that the parallel code solves the same
problem as the serial one. The validation, which we call experiment~0,
consists of twelve small cases. They differ in grid size (from $3\times3$ to
$40\times32$ cells), in the number of steps (8 to 50), in the layout (1 to 6
ranks with 1 to 4 threads each), in the decomposition (row-wise or
two-dimensional) and in the halo exchange (blocking or overlapped). Ten cases
carry a pigment drop in the wind-driven flow, one of them under a periodic
wind. The other two have a known answer: a constant field, and a cosine mode
that decays in a fluid at rest.

Each case is run by three programs: the serial solver \texttt{heat\_serial},
the same solver compiled with OpenMP and run on one thread
(\texttt{heat\_serial\_omp}), and the parallel solver \texttt{heat\_parallel}.
We keep the second serial reference because the speedups of the later
experiments are measured against a serial program, and the two are only
comparable if they solve the same problem. The following checks are applied.
\begin{enumerate}
  \item \textbf{Mass conservation} (every case). The scheme conserves the
  total pigment $M=\Delta x\,\Delta y\sum_{i,j}c_{ij}$ exactly, so the relative
  change of $M$ over the serial run must stay below $10^{-12}$.
  \item \textbf{Serial references} (every case). The final field of
  \texttt{heat\_serial\_omp} is compared with that of \texttt{heat\_serial},
  cell by cell.
  \item \textbf{Serial against parallel} (every case). The final field of
  \texttt{heat\_parallel} is compared with the serial one, cell by cell. A
  case passes only if the difference is at most $10^{-10}+10^{-10}\,|c|$ in
  every cell; the second check uses the same tolerance.
  \item \textbf{Constant field} (one case). A uniform concentration must stay
  uniform in the moving fluid, to within $10^{-12}$ of its value.
  \item \textbf{Cosine decay} (one case). The numerical result is compared
  with the exact solution of the discrete scheme, with a tolerance of
  $10^{-12}$. The error against the continuous solution is the discretization
  error, which we record but do not use as a pass criterion.
\end{enumerate}
This gives 38 checks with a pass criterion. Their results are reported in
Section~\ref{sec:res-correctness}.
```

### 4.2 Flow and Pigment Evolution (replace the whole subsection)

The old text described the picture without saying why it is in a validation section. The new text ties it to the equations of Section 2: each thing we see is something the model should do. The figure is new, on a $128\times128$ grid instead of $64\times64$, with the same parameters and the same times $\tau$.

```latex
\subsection{Flow and Pigment Evolution}

The checks of Section~\ref{sec:validation} compare numbers on small grids. As
a complement, we look at one full run and ask whether the solution behaves as
the equations of Section~2 say it should. The run uses the unit square with a
$128\times128$ grid, the steady wind of amplitude $F_0=0.8$, viscosity
$\nu=0.004$ and pigment diffusivity $\kappa=10^{-4}$. The flow starts from
rest and spins up for 8 time units before the pigment blob of Section~2.3 is
released; we write $\tau$ for the time since the release.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{flow_evolution}
\caption{Wind-driven circulation and pigment transport in the unit square on
a $128\times128$ grid, with $F_0=0.8$, $\nu=0.004$ and $\kappa=10^{-4}$.
(a) Streamlines at the release of the pigment, after a spin-up of 8 time
units; the dashed circle marks the blob. (b)--(d) Pigment concentration at
times $\tau$ after the release, on a common square-root colour scale, with
the streamlines in the background.}
\label{fig:flow}
\end{figure}

Figure~\ref{fig:flow}(a) shows the flow at the release. The forcing
$F_0\sin(2\pi y/L_y)$ changes sign at mid-height, so it drives two gyres that
turn in opposite directions, one in each half of the basin, and the
streamlines close inside the basin because the walls are impermeable. Between
the gyres both of them move the fluid to the left, and this is where the blob
is placed.

The pigment then follows the flow. It is first carried to the left wall,
where it splits into two branches, one for each gyre. By $\tau=2.44$ the
branches have turned along the top and bottom walls
(Figure~\ref{fig:flow}(c)), and by $\tau=9.77$ the pigment has travelled
around both gyres and diffuses towards their centres
(Figure~\ref{fig:flow}(d)). Meanwhile the peak concentration falls from 1 to
0.26 and then to 0.06, while the total pigment mass stays the same to within
rounding, as required by Equation~(12). The run thus shows advection along
the streamlines, diffusion across them and no loss through the walls, which
is the behaviour the model is meant to have.
```

---

## Section 5

Changes: shorter throughout, 5.3.3 split into two subsections, order aligned with Section 6, and the setup that sat in Section 6 moved here. Tables 6 and 7 (weak-scaling grids, GPU grids) are replaced by one sentence each, because Tables 12 and 18 already list the same grids.

### 5.1 Computing Environment (replace)

```latex
\subsection{Computing Environment}
\label{sec:environment}

All experiments ran on the NUS Vanda cluster (Table~\ref{tab:environment}).
The CPU runs use a Release build with GCC and Open MPI, and the GPU runs an
NVHPC build with OpenMP target offloading of the CG solver. Each experiment
is configured by a JSON file and submitted as a PBS job that reserves whole
nodes.

% unchanged: Table 2, "Computing environment for the performance experiments"

Every thread is bound to a physical core of its own. Unless stated otherwise,
a run of at most 32 cores stays on socket~0, and a larger run is split evenly
between the two sockets of each node, with at most 32 cores per socket and 64
per node. Runs of 128 and 256 cores therefore span two and four nodes.
```

### 5.2 Workloads and Performance Measurement (replace)

```latex
\subsection{Workloads and Performance Measurement}
\label{sec:measurement}

The fixed-size experiments share one benchmark problem, which we call $P$: a
$1024\times1024$ grid on the unit square, run for 100 steps from rest
(Table~\ref{tab:benchmark}).

% unchanged: Table 3, "Fixed-size benchmark parameters and measurements"

\paragraph{Timings.} We time the whole time-stepping loop and write $T$ for
its wall-clock time. In a parallel run every rank times its own loop and $T$
is the maximum over the ranks. Spin-up and output lie outside the loop.
Within the loop we record the time $T_P$ spent in the Poisson solves and the
mean number $\bar k$ of CG iterations per step, which give the time of one CG
iteration,
\begin{equation}
  t_{\mathrm{it}} = \frac{T_P}{N_s\,\bar k},
\end{equation}
with $N_s$ the number of steps. The remainder $T-T_P$ is the pigment update,
the vorticity update and the velocity calculation. Inside the Poisson solves
we also record the time in the halo exchanges and in the global sums
(\texttt{MPI\_Allreduce}), and report both per CG iteration as \emph{halo}
and \emph{sum}. These two are maxima over the ranks, so they do not add up to
$t_{\mathrm{it}}$. The follow-up runs also record the minimum of the
global-sum time over the ranks (\emph{min sum}): the time of the slowest rank
includes its waiting for the others, while that of the fastest rank is close
to the cost of the reduction alone.

\paragraph{Repetitions.} Unless stated otherwise, a configuration has one
warm-up run followed by five measured repetitions. We report the median and,
as the spread, the largest minus the smallest value.

\paragraph{Speedup and efficiency.} For a run on $p$ cores,
\begin{equation}
  S(p) = \frac{T_s}{T(p)}, \qquad E(p) = \frac{S(p)}{p},
  \label{eq:speedup}
\end{equation}
where $T_s$ is the loop time of the serial program \texttt{heat\_serial} on
one core. The second serial reference \texttt{heat\_serial\_omp}
(Section~\ref{sec:validation}) is timed alongside it.
```

### 5.3 Experimental Design (replace, including all subsections up to GPU Offloading)

```latex
\subsection{Experimental Design}

The experiments are listed in the order of their results in
Section~\ref{sec:results}. All of them use the measurements of
Section~\ref{sec:measurement}, and all parallel runs use the two-dimensional
decomposition with blocking halo exchange unless stated otherwise.

\subsubsection{Serial Runtime and Problem Size}
\label{sec:setup-serial}

\texttt{heat\_serial} runs on one core for grids of $256^2$, $512^2$,
$1024^2$, $2048^2$ and $4096^2$ cells on the unit square, with 100 steps and
the same physical parameters and CG tolerance. The time step follows the
stability limit of each grid. The $256^2$ and $512^2$ grids have five
repetitions, $1024^2$ and $2048^2$ three, and $4096^2$ one. An $n\times n$
grid has $N=(n-1)^2$ streamfunction unknowns, and the five-point operator,
which is applied without storing a matrix, corresponds to
\begin{equation}
  \mathrm{nnz} = 5N - 4(n-1)
\end{equation}
nonzero entries. We report the time per unknown and iteration,
$t_{\mathrm{it}}/N$.

\subsubsection{MPI Strong Scaling}
\label{sec:setup-strong-mpi}

Experiment~1 solves $P$ on
\begin{equation}
  p \in \{1, 2, 4, 8, 16, 32, 64, 128, 256\}
\end{equation}
ranks with one thread each. Up to 64 ranks fit on one node; 128 and 256 ranks
use two and four nodes.

\subsubsection{OpenMP Strong Scaling}
\label{sec:setup-openmp}

Experiment~2 solves $P$ on one rank with
\begin{equation}
  p \in \{1, 2, 4, 8, 16, 32, 36, 64, 72\}
\end{equation}
threads on one node. The first 36 threads fill socket~0 and the rest go to
socket~1, so 36 threads use one whole socket and 72 the whole node. We stop
at 72 because more threads would put two of them on one core. The runs with 2
to 72 threads were repeated after we introduced parallel first-touch
initialization, in which every thread zeroes its own part of the arrays so
that the memory is placed on its socket. We report these later runs.

\subsubsection{MPI Weak Scaling}
\label{sec:setup-weak}

Experiment~5 gives every rank a block of $256\times256$ cells and one thread.
With $p=q^2$ ranks on a $q\times q$ process grid, the global grid has
$256q\times256q$ cells on a $q\times q$ basin, so the grid spacing stays at
$1/256$; the physical parameters, the time-step settings and the 100 steps
are those of $P$. We use
\begin{equation}
  p \in \{1, 4, 9, 16, 36, 64, 144, 256\},
\end{equation}
all in one job on the same four nodes. Up to 64 ranks run on one node, 144 on
three and 256 on four. Up to 16 ranks run on socket~0 only; from 36 ranks on,
the ranks are split evenly between the two sockets of each node. We measure
two efficiencies against the one-rank run,
\begin{equation}
  E_w(p) = \frac{T(1)}{T(p)}, \qquad
  E_{w,\mathrm{it}}(p) = \frac{t_{\mathrm{it}}(1)}{t_{\mathrm{it}}(p)}.
\end{equation}

A placement control on one node separates what changes between 16 and 64
ranks. It keeps the $256\times256$ block per rank and runs 16 ranks on one
socket or as $8+8$, 32 ranks on one socket or as $16+16$, and 64 ranks as
$32+32$ twice, the second time with the clock of all 72 cores logged every
5~s.

\subsubsection{Hybrid MPI--OpenMP Layouts}
\label{sec:setup-hybrid}

Experiment~3 solves $P$ on 256 cores of four nodes, 64 per node, and varies
how they are split between ranks and threads (Table~\ref{tab:layouts}). These
runs have ten repetitions and use first-touch initialization.

% unchanged: Table 4, "MPI--OpenMP layouts for the fixed 256-core comparison"

Two controls make the layouts comparable. The first splits the 64 cores of
one node between ranks and threads, from $64\times1$ to $2\times32$, so that
no message leaves the node. The second keeps 64 ranks with one thread each
and spreads them over one, two and four nodes, for $P$ and for the
$2048\times2048$ problem of the weak scaling.

\subsubsection{Decomposition and Halo Exchange}
\label{sec:setup-comm}

Experiment~4 solves $P$ on 256 ranks with one thread each, on four nodes,
with the four combinations of Table~\ref{tab:comm-configs}. The row-wise
decomposition gives each rank a strip of 4 rows of 1024 cells; the
two-dimensional one gives it a $64\times64$ block on a $16\times16$ process
grid. The chosen exchange is used for both the pigment update and the CG
matrix--vector product.

% unchanged: Table 5, "Decomposition and halo exchange configurations using 256 CPU cores"

\subsubsection{GPU Offloading}
\label{sec:setup-gpu}

In the GPU runs only the CG loops run on the A40. The pigment, vorticity and
velocity updates stay on the CPU threads of the rank, and the global sums and
halo exchanges go through the host. We test two layouts: one rank of 36
threads (one socket) with one A40, and two ranks of 36 threads (the whole
node) with two A40, one GPU per rank. Each is compared with the CPU run of
the same layout, on grids of $512^2$ and $1024^2$ cells for 100 steps,
$2048^2$ for 20 steps and $4096^2$ for 10 steps. Because the step counts
differ, we compare the time per CG iteration, and we also record the
host--device copy time. The CPU and GPU runs use different compilers and MPI
libraries (Table~\ref{tab:environment}).

% unchanged: "Three-Dimensional Experiments", now 5.3.8
```

---

## Section 6

### Opening of Section 6 (replace the three paragraphs before 6.0)

The JSON and PBS sentence moved to 5.1, and the list of experiments is already in 5.3.

```latex
\section{Results and Performance Analysis}
\label{sec:results}

As the central piece of our concerned topic, this section presents what we
find from the experiments of Section~5. We start from the correctness check
and the serial runtime, proceed to strong scaling with MPI and with OpenMP
and to weak scaling with MPI, and then turn to the hybrid layouts and the
communication strategies, before closing with the GPU and the
three-dimensional extension.

Apart from the planned runs, several follow-up runs appear along the way.
They were motivated by some unexpected behaviors in the first batch of
results. For example, in the weak scaling, the time per conjugate gradient
(CG) iteration jumps as the number of ranks grows, and the jumps could have
different roots of cause, such as a lower CPU clock, the density of ranks
within a socket, or the cost of communication between nodes. The extra runs
separate these factors.
```

### 6.0 Correctness Check (replace the whole subsection)

```latex
\subsection{Correctness Check}
\label{sec:res-correctness}

All 38 checks of Section~\ref{sec:validation} pass
(Table~\ref{tab:validation}). The mass changes of the serial runs are at most
$8\times10^{-16}$, which is far below the tolerance. The fields of the two
serial programs agree exactly, with maximum differences found to be zero
across the 12 cases, so either of them can serve as the reference. Between
the serial and the parallel fields, the largest difference found is
$9.5\times10^{-13}$, with nine cases ranging between $1\times10^{-14}$ and
$1\times10^{-12}$, and the rest holding zero, namely a full precision. We
conclude that these nonzero values are rounding differences, as the global
sums of the conjugate gradient method add the same terms in a different order
across ranks and threads.

\begin{table}[htbp]
\centering
\small
\begin{tabular}{llrr}
\toprule
Check & Measured quantity & Tolerance & Largest value \\
\midrule
Mass conservation & relative mass change & $10^{-12}$ & $8.04\times10^{-16}$ \\
Serial references & field difference & $10^{-10}$ & $0$ \\
Serial against parallel & field difference & $10^{-10}$ & $9.53\times10^{-13}$ \\
Constant field & deviation from the initial value & $10^{-12}$ (relative) & $1.42\times10^{-14}$ \\
Cosine decay & error against the discrete solution & $10^{-12}$ & $9.99\times10^{-16}$ \\
\bottomrule
\end{tabular}
\caption{Results of the correctness check over the twelve cases. The error of
the cosine case against the continuous solution, $7.89\times10^{-8}$, is the
discretization error and has no pass criterion.}
\label{tab:validation}
\end{table}
```

### 6.1 Problem Size versus Serial Runtime

No change. It has no setup that repeats Section 5.

### 6.2 MPI Strong Scaling (replace the two paragraphs before Table 10)

The benchmark, the repetitions and the definitions of $T$, $T_P$, halo, sum and $t_{\mathrm{it}}$ moved to 5.2. Everything from Table 10 on stays as it is.

```latex
\subsection{MPI Strong Scaling}

Now, we may proceed to discuss how strong scaling works with MPI on the
benchmark problem by varying the rank count. The reference $T_s$ is the pure
serial program \texttt{heat\_serial} (273.7~s; \texttt{heat\_serial\_omp}
takes 274.0~s, within 0.2\%). In this experiment, we concerned most about the
cost of one iteration $t_{\mathrm{it}}$, since it separates the parallel
implementation from the algorithm, which demands iterations on a bigger grid
(see our discussion in Section~\ref{sec:res-serial}). Apart from the Poisson
solves, the remaining parts together take an unassuming amount of time
$T-T_P$, and a corresponding analysis of it will be needed only if it ever
grows to a sizable part of $T$.
```

### 6.3 OpenMP Strong Scaling (replace the two paragraphs before "We will elaborate every potential source")

The sentences on the 72-thread limit and the binding protocol moved to 5.3.3.

```latex
\subsection{OpenMP Strong Scaling}

We now study the strong scaling with OpenMP, where multiple threads of a
single rank work on the same benchmark problem together. Because the threads
share memory, we expect different dynamics from strong scaling with MPI.
There, the gain from splitting the domain among ranks competes with the
growing cost of communication. Here, the limit comes from the cost of the
threads working together, namely the reductions and the barriers between
them.

According to Figure~\ref{fig:strong-openmp}, the runtime drops as the number
of threads increases, close to the ideal estimate. From 32 threads on, it
stops improving and both runtime and the speedup start to deviate, yet it did
not worsen dramatically. So, the cost that halts the scaling does not grow
with the number of threads and acts like a floor. This is confirmed with a
stable count of $t_{\mathrm{it}}$ between 32 and 72 threads.
```

### 6.4 MPI Weak Scaling (replace the whole subsection)

The setup paragraphs moved to 5.3.4. Table 12 loses the spread, nodes and halo columns. The new figure and its paragraph come before the one-node and multi-node discussion, because the figure is what splits the analysis into those two parts.

```latex
\subsection{MPI Weak Scaling}

In this section, we start exploring how a weak scaling over MPI ranks
behaves. It is worth noting that since we need to enlarge the problem in
order to assign the same $256^2$ block to each rank, there exists a nonlinear
growth in the required number of iterations. This as a result makes the cost
per iteration $t_{\mathrm{it}}$ a proper measure of efficiency compared to
the total runtime $T$. Hence, on top of the conventional efficiency $E_w$, we
analyze $E_{w,\mathrm{it}}$ and its dynamics.

% Table 12 with nine columns (same as results/figures/weak_scaling_mpi_table.tex)
\begin{table}[!t]
\centering
\small
\begin{tabular}{rrrrrrrrr}
\toprule
Ranks $p$ & Grid & $T$ (s) & $E_w$ & $\bar k$ & $t_{\mathrm{it}}$ (ms) & $E_{w,\mathrm{it}}$ & sum (ms) & min sum (ms) \\
\midrule
serial & $256^2$ & 3.88 & -- & -- & -- & -- & -- & -- \\
\midrule
1 & $256^2$ & 3.91 & 1.000 & 337 & 0.114 & 1.00 & 0.000 & 0.000 \\
4 & $512^2$ & 8.65 & 0.453 & 687 & 0.125 & 0.91 & 0.008 & 0.005 \\
9 & $768^2$ & 13.34 & 0.293 & 1038 & 0.128 & 0.90 & 0.012 & 0.006 \\
16 & $1024^2$ & 18.32 & 0.214 & 1390 & 0.131 & 0.87 & 0.014 & 0.009 \\
36 & $1536^2$ & 29.33 & 0.133 & 2095 & 0.140 & 0.82 & 0.021 & 0.009 \\
64 & $2048^2$ & 47.77 & 0.082 & 2801 & 0.170 & 0.67 & 0.026 & 0.012 \\
144 & $3072^2$ & 133.42 & 0.029 & 4216 & 0.316 & 0.36 & 0.182 & 0.123 \\
256 & $4096^2$ & 228.30 & 0.017 & 5634 & 0.405 & 0.28 & 0.261 & 0.184 \\
\bottomrule
\end{tabular}
\caption{Experiment 5: weak scaling over MPI ranks, medians of five repetitions. $E_w=T(1)/T(p)$ and $E_{w,\mathrm{it}}=t_{\mathrm{it}}(1)/t_{\mathrm{it}}(p)$ are taken against the one-rank run. Sum is the time per CG iteration inside the global sums for the slowest rank and min sum for the fastest rank. Up to 64 ranks the runs are on one node; 144 ranks use three nodes and 256 ranks four.}
\label{tab:weak-mpi}
\end{table}

According to Table~\ref{tab:weak-mpi}, $E_w$ falls to 0.45 at four ranks and
to 0.017 at 256, but most of this loss is algorithmic rather than parallel.
Since the Poisson solves take nearly the whole loop,
$T(p)\approx N_s\,\bar k(p)\,t_{\mathrm{it}}(p)$, and therefore
\begin{equation}
  E_w(p) = \frac{T(1)}{T(p)} \approx
  \frac{\bar k(1)}{\bar k(p)}\cdot\frac{t_{\mathrm{it}}(1)}{t_{\mathrm{it}}(p)}
  = \frac{\bar k(1)}{\bar k(p)}\,E_{w,\mathrm{it}}(p).
  \label{eq:weak-factor}
\end{equation}
The first factor is set by the algorithm. The grid has $n=256\sqrt p$ cells a
side, and the number of iterations of unpreconditioned CG grows in proportion
to $n$: $\bar k/\sqrt p$ stays between 337 and 352 for all eight points, and
at 256 ranks $\bar k$ is 16.7 times its one-rank value. The second factor is
the parallel efficiency proper. At 256 ranks the two give
$0.28/16.7=0.017$, the measured $E_w$.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{weak_scaling_mpi_breakdown_summin}
\caption{Experiment 5: breakdown of the weak-scaled runs. Top: the Poisson
solves in the loop time, as a share and in seconds. Bottom: one CG iteration
split into the halo exchange, the global sums and the rest, which is the
local computation, as a share and in milliseconds. The global sums are those
of the fastest rank.}
\label{fig:weak-breakdown}
\end{figure}

Figure~\ref{fig:weak-breakdown} shows where the time goes. The Poisson solves
take 98.6\% to 99.9\% of the loop at every rank count, so again it is the
iteration that we need to explain. We split it into the halo exchange, the
global sums and the rest, namely the local computation. The sums are taken
from the fastest rank, since the time of the slowest one also contains its
waiting for the others. Two regimes are then clear. Within one node, the
computation holds at least 87\% of an iteration, and it is the computation
itself that grows, from 0.114~ms at one rank to 0.148~ms at 64. Beyond one
node, the computation stays at 0.14 to 0.15~ms, yet the halo exchange and the
global sums now take 55\% of the iteration at 144 ranks and 64\% at 256. So
the efficiency lost within a node and the efficiency lost across nodes have
different origins, and we discuss them in turn.

Within one node, $E_{w,\mathrm{it}}$ declines gently to 0.87 at 16 ranks and
0.82 at 36, then drops to 0.67 at 64. The placement matters here: up to 16
ranks run on socket~0 only, while 36 and 64 ranks are split evenly between
the two sockets, with 18 and 32 ranks per socket. The control runs of
Table~\ref{tab:placement} separate the possible causes. Moving 16 ranks from
one socket to $8+8$ leaves $t_{\mathrm{it}}$ unchanged (0.132~ms in both
cases), so crossing between the sockets costs nothing. Putting 32 ranks on
one socket instead of $16+16$ raises $t_{\mathrm{it}}$ from 0.142 to
0.173~ms, close to the 0.180~ms of 64 ranks with 32 per socket. Within a
node, the time per iteration thus follows the number of ranks per socket, not
the total number of ranks. The CPU clock stayed at 2.5~GHz on the busy cores
during the second 64-rank control run (median 2500~MHz, 5th percentile
2438~MHz), which rules out thermal throttling. The ranks of one socket share
its last-level cache and its memory channels. Each rank's CG works on five
arrays of about 0.53~MB, so 16 ranks need about 43~MB and 32 ranks about
85~MB, against an L3 cache of 67.5~MB per socket. A combined working set that
no longer fits in the cache is consistent with the step at 32 ranks per
socket, but we did not measure cache misses or memory traffic, so this
explanation remains a hypothesis.

% unchanged: Table 13, "Placement control on one node"

Beyond one node, $E_{w,\mathrm{it}}$ falls to 0.36 on three nodes (144 ranks)
and 0.28 on four (256 ranks), although these runs place fewer ranks per
socket (24 at 144 ranks) or as many (32 at 256 ranks) as the one-node run
with 64 ranks. The global sums account for the difference: at 256 ranks the
slowest rank spends 0.261~ms per iteration in them and even the fastest
spends 0.184~ms, against 0.026 and 0.012~ms with 64 ranks on one node. So the
reduction between the nodes, not the waiting for slower ranks, dominates the
iteration. The node-spread control of Section~\ref{sec:res-hybrid}
(Table~\ref{tab:node-spread}) shows the same effect with the rank count fixed
at 64.

All points of Table~\ref{tab:weak-mpi} come from one job on the same four
nodes. An earlier run of the 144-rank point on other nodes took 156.8~s
instead of 133.4~s, so multi-node times differ by about 15\% between jobs,
more than the spread within one job (5\% at 256 ranks).
```

### 6.5 Hybrid Layouts (replace the first paragraph only)

I suggest the title "Hybrid MPI--OpenMP Layouts", since OpenMP strong scaling is 6.3. The layouts, the four nodes and the two controls are now described in 5.3.5. The rest of the subsection stays as it is.

```latex
\subsection{Hybrid MPI--OpenMP Layouts}
\label{sec:res-hybrid}

In this section, we aim to test how results would be different by changing
how the number of cores split between ranks and threads, with the total
fixed at 256 cores on four Vanda nodes. Realized that these seven layouts is
only informative up to a sense of speed and any further comparisons is hard
unless comparable setups are provided, we also draw on the two control tests
of Section~\ref{sec:setup-hybrid}. These extras can then answer questions
such as how varying the layout within a node would produce.
```

### 6.6 Decomposition and Communication Overlap (replace the first paragraph only)

The geometry of the strips and blocks moved to 5.3.6.

```latex
\subsection{Decomposition and Communication Overlap}

Table~\ref{tab:comm} shows the results of the four combinations of
decomposition and halo exchange on 256 ranks
(Section~\ref{sec:setup-comm}).
```

### 6.7 GPU-Assisted Whole-Simulation Performance (replace the first paragraph only)

What runs on the GPU, the layouts and the compilers moved to 5.3.7. The two sentences that stay are checks on the measurements, not setup.

```latex
\subsection{GPU-Assisted Whole-Simulation Performance}

Table~\ref{tab:gpu} compares each GPU configuration with the CPU run of the
same layout (Section~\ref{sec:setup-gpu}). The node type does not change CPU
timings: the CPU code with one rank of 36 threads took 11.24~s on a GPU node
(experiment~2) and 11.19~s on a CPU node. The GPU timings are reproducible:
the $1024^2$ runs of 3 and 8 November differ by less than 1\% (31.80 against
31.49~s on one GPU).
```

---

## Noticed but not changed

- **Table 9 has no $N$ or nnz column.** Section 5 promises the number of unknowns and nonzeros for the serial sizes, but the table does not show them. The values are below if you want to add two columns.

  | Grid | $N=(n-1)^2$ | nnz $=5N-4(n-1)$ |
  | --- | --- | --- |
  | $256^2$ | 65 025 | 324 105 |
  | $512^2$ | 261 121 | 1 303 561 |
  | $1024^2$ | 1 046 529 | 5 228 553 |
  | $2048^2$ | 4 190 209 | 20 942 857 |
  | $4096^2$ | 16 769 025 | 83 828 745 |

- **Section 4.2 has two hard-coded references**, "Section 2.3" and "Equation (12)". Replace them with `\ref` to your labels.
- **Section 5.3.8 (3D)** still says "experiments 1, 4 and 5" and "problem $P$", which now match the names used in 5.2 and 5.3.
- **Sections 6.6 to 6.8** keep their analysis text. Only their opening setup paragraphs changed.
