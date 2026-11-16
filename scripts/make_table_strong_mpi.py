"""LaTeX table of experiment 1 (strong scaling over MPI ranks) from results/raw."""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "results/raw"
groups = defaultdict(list)
for name in ("strong_scaling_mpi_1node.csv", "strong_scaling_mpi_4nodes.csv"):
    for r in csv.DictReader((RAW / name).open(newline="")):
        if r["status"] == "ok" and r["phase"] == "measured":
            groups[r["run_id"]].append(r)


def med(rows, key):
    return st.median(float(r[key]) for r in rows)


serial = med(groups["serial_reference"], "elapsed_s")
lines = []
t = [float(r["elapsed_s"]) for r in groups["serial_reference"]]
lines.append(f"serial & {st.median(t):.1f} & {max(t) - min(t):.1f} & " + " & ".join(["--"] * 7) + r" \\")
lines.append(r"\midrule")
for p in (1, 2, 4, 8, 16, 32, 64, 128, 256):
    rows = groups[f"p{p}_t1_2d_blocking"]
    t = [float(r["elapsed_s"]) for r in rows]
    T = st.median(t)
    steps = float(rows[0]["steps"])
    k = med(rows, "poisson_iterations_mean")
    n = steps * k
    tp, halo, tsum = (med(rows, c) for c in ("poisson_s", "poisson_halo_s", "poisson_sum_s"))
    S = serial / T
    lines.append(f"{p} & {T:.2f} & {max(t) - min(t):.2f} & {S:.1f} & {S / p:.2f} & "
                 f"{100 * tp / T:.1f} & {k:.0f} & {1e3 * tp / n:.3f} & {1e3 * halo / n:.3f} & {1e3 * tsum / n:.3f} \\\\")
head = r"""\begin{table}[htbp]
\centering
\small
\begin{tabular}{rrrrrrrrrr}
\toprule
Ranks $p$ & $T$ (s) & spread (s) & $S$ & $E$ & $T_P/T$ (\%) & $\bar k$ & $t_{\mathrm{it}}$ (ms) & halo (ms) & sum (ms) \\
\midrule
"""
tail = r"""\bottomrule
\end{tabular}
\caption{Experiment 1: strong scaling over MPI ranks on the fixed problem $P$ ($1024\times1024$ cells, 100 steps), one thread per rank. $T$ is the median of five repetitions of the loop time and the spread is their largest minus their smallest value; $S=T_{\mathrm{serial}}/T$ and $E=S/p$ are taken against \texttt{heat\_serial}. $T_P/T$ is the share of the Poisson solves, $\bar k$ the mean CG iterations per step, and $t_{\mathrm{it}}$, halo and sum the time per CG iteration of the Poisson solve, of its halo exchanges and of its global sums. Runs up to 64 ranks are on one node; 128 and 256 ranks use two and four nodes.}
\label{tab:strong-mpi}
\end{table}
"""
out = head + "\n".join(lines) + "\n" + tail
(Path(__file__).resolve().parents[1] / "results/figures/strong_scaling_mpi_table.tex").write_text(out)
print(out)
