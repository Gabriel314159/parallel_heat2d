"""LaTeX tables of the two controls of experiment 3 (node spread, one-node split) from results/raw."""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results/raw"


def load(name):
    return [r for r in csv.DictReader((RAW / name).open(newline="")) if r["status"] == "ok" and r["phase"] == "measured"]


def med(rows, key):
    return st.median(float(r[key]) for r in rows)


def per_it(rows, key):
    """Median of a time column in ms per CG iteration (key None: t_it - sum)."""
    n = float(rows[0]["steps"]) * med(rows, "poisson_iterations_mean")
    if key is None:
        return 1e3 * (med(rows, "poisson_s") - med(rows, "poisson_sum_s")) / n
    return 1e3 * med(rows, key) / n


REST = r"\shortstack[r]{$t_{\mathrm{it}}-\mathrm{sum}$\\(ms)}"

# ---- 64 ranks on 1, 2 and 4 nodes ---------------------------------------------------------------------------
groups = defaultdict(list)
for r in load("node_spread.csv"):
    groups[(int(r["nx"]), int(r["nodes"]))].append(r)
lines = []
for nx in sorted({k[0] for k in groups}):
    if lines:
        lines.append(r"\midrule")
    base = med(groups[(nx, 1)], "elapsed_s")
    for nodes in sorted(k[1] for k in groups if k[0] == nx):
        rows = groups[(nx, nodes)]
        T = med(rows, "elapsed_s")
        lines.append(f"${nx}^2$ & {nodes} & {T:.2f} & {T / base:.2f} & {per_it(rows, 'poisson_sum_s'):.3f} & "
                     f"{per_it(rows, None):.3f} \\\\")
spread = r"""\begin{table}[htbp]
\centering
\small
\begin{tabular}{rrrrrr}
\toprule
Problem & Nodes & \shortstack[r]{$T$\\(s)} & $\dfrac{T}{T_{1\,\mathrm{node}}}$ & \shortstack[r]{sum\\(ms)} & """ + REST + r""" \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\caption{The same 64 MPI ranks, one thread each, on one, two and four nodes (64, 32 and 16 ranks per node, split evenly between the two sockets); 2D decomposition, blocking exchange, 100 steps. $1024^2$ is the fixed problem $P$; $2048^2$ is the 64-rank problem of experiment 5, a $256\times256$ block per rank. $T$ is the median of five repetitions of the loop time and $T/T_{1\,\mathrm{node}}$ its ratio to the same problem on one node. sum is the time per CG iteration of the global sums (maximum over the ranks) and $t_{\mathrm{it}}-\mathrm{sum}$ the rest of the iteration: the local sweeps and the halo exchange.}
\label{tab:node-spread}
\end{table}
"""

# ---- the 64 cores of one node -----------------------------------------------------------------------------
groups = defaultdict(list)
for r in load("hybrid_node.csv"):
    groups[(int(r["mpi_ranks"]), int(r["omp_threads"]))].append(r)
layouts = sorted(groups, key=lambda k: -k[0])
base = med(groups[layouts[0]], "elapsed_s")  # 64 x 1
lines = []
for p, t in layouts:
    rows = groups[(p, t)]
    T = med(rows, "elapsed_s")
    lines.append(f"${p}\\times{t}$ & {T:.2f} & {base / T:.2f} & {per_it(rows, None):.3f} \\\\")
node = r"""\begin{table}[htbp]
\centering
\small
\begin{tabular}{rrrr}
\toprule
\shortstack[r]{Ranks $\times$\\threads} & \shortstack[r]{$T$\\(s)} & $\dfrac{T_{\mathrm{MPI}}}{T}$ & """ + REST + r""" \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\caption{The 64 cores of one node (32 per socket) split between MPI ranks and OpenMP threads, on the fixed problem $P$ ($1024\times1024$ cells, 100 steps), 2D decomposition, blocking exchange; no message leaves the node. $T$ is the median of five repetitions of the loop time and $T_{\mathrm{MPI}}/T$ the speedup over the pure MPI layout $64\times1$. $t_{\mathrm{it}}-\mathrm{sum}$ is the time per CG iteration without its global sums: the local sweeps and the halo exchange.}
\label{tab:hybrid-node}
\end{table}
"""

(ROOT / "results/figures/node_spread_table.tex").write_text(spread)
(ROOT / "results/figures/hybrid_node_table.tex").write_text(node)
print(spread)
print(node)
