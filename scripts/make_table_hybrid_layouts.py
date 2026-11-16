"""LaTeX table of experiment 3 (pure MPI vs MPI + OpenMP, 256 cores) from results/raw."""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "results/raw"
groups = defaultdict(list)
for r in csv.DictReader((RAW / "hybrid_layouts_firsttouch.csv").open(newline="")):
    if r["status"] == "ok" and r["phase"] == "measured":
        groups[(int(r["mpi_ranks"]), int(r["omp_threads"]))].append(r)


def med(rows, key):
    return st.median(float(r[key]) for r in rows)


layouts = sorted(groups, key=lambda k: -k[0])
base = med(groups[layouts[0]], "elapsed_s")  # 256 x 1, pure MPI
lines = []
for p, t in layouts:
    rows = groups[(p, t)]
    r0 = rows[0]
    T = med(rows, "elapsed_s")
    cells = f"{int(r0['nx']) // int(r0['px'])}\\times{int(r0['ny']) // int(r0['py'])}"
    k = med(rows, "poisson_iterations_mean")
    n = float(r0["steps"]) * k
    tp, halo, tsum = (med(rows, c) for c in ("poisson_s", "poisson_halo_s", "poisson_sum_s"))
    lines.append(f"${p}\\times{t}$ & ${cells}$ & {T:.2f} & {base / T:.2f} & "
                 f"{1e3 * tp / n:.3f} & {1e3 * halo / n:.3f} & {1e3 * tsum / n:.3f} & {1e3 * (tp - tsum) / n:.3f} \\\\")

out = r"""\begin{table}[htbp]
\centering
\small
\begin{tabular}{rrrrrrrr}
\toprule
\shortstack[r]{Ranks $\times$\\threads} & \shortstack[r]{Cells\\per rank} & \shortstack[r]{$T$\\(s)} & $\dfrac{T_{\mathrm{MPI}}}{T}$ & \shortstack[r]{$t_{\mathrm{it}}$\\(ms)} & \shortstack[r]{halo\\(ms)} & \shortstack[r]{sum\\(ms)} & \shortstack[r]{$t_{\mathrm{it}}-\mathrm{sum}$\\(ms)} \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\caption{Experiment 3: the 256 cores of four nodes split between MPI ranks and OpenMP threads, on the fixed problem $P$ ($1024\times1024$ cells, 100 steps), 2D decomposition, blocking exchange. Cells per rank are $n_x\times n_y$. $T$ is the median of ten repetitions of the loop time and $T_{\mathrm{MPI}}/T$ the speedup over the pure MPI layout $256\times1$. The Poisson solves are 99.7 to 100\% of $T$ and take $\bar k=1095$ CG iterations per step in every layout; $t_{\mathrm{it}}$, halo and sum are the time per CG iteration of the solve, of its halo exchanges and of its global sums (maxima over the ranks, so the parts do not add up); $t_{\mathrm{it}}-\mathrm{sum}$ is the rest of the iteration, the local sweeps and the halo exchange.}
\label{tab:hybrid-layouts}
\end{table}
"""
(Path(__file__).resolve().parents[1] / "results/figures/hybrid_layouts_table.tex").write_text(out)
print(out)
