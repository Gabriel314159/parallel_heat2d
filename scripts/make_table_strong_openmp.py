"""LaTeX table of experiment 2 (strong scaling over OpenMP threads) from results/raw."""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "results/raw"


def load(name):
    groups = defaultdict(list)
    for r in csv.DictReader((RAW / name).open(newline="")):
        if r["status"] == "ok" and r["phase"] == "measured":
            groups[r["run_id"]].append(r)
    return groups


def med(rows, key):
    return st.median(float(r[key]) for r in rows)


ref, runs = load("strong_scaling_openmp.csv"), load("strong_scaling_openmp_firsttouch.csv")
t = [float(r["elapsed_s"]) for r in ref["serial_reference"]]
serial = st.median(t)
lines = [f"serial & {serial:.1f} & {max(t) - min(t):.1f} & " + " & ".join(["--"] * 7) + r" \\", r"\midrule"]
for p in (2, 4, 8, 16, 32, 36, 64, 72):
    rows = runs[f"p1_t{p}_2d_blocking"]
    t = [float(r["elapsed_s"]) for r in rows]
    T = st.median(t)
    k = med(rows, "poisson_iterations_mean")
    n = float(rows[0]["steps"]) * k
    tp, halo, tsum = (med(rows, c) for c in ("poisson_s", "poisson_halo_s", "poisson_sum_s"))
    S = serial / T
    lines.append(f"{p} & {T:.2f} & {max(t) - min(t):.2f} & {S:.1f} & {S / p:.2f} & "
                 f"{100 * tp / T:.1f} & {k:.0f} & {1e3 * tp / n:.3f} & {1e3 * halo / n:.3f} & {1e3 * tsum / n:.4f} \\\\")
head = r"""\begin{table}[htbp]
\centering
\small
\begin{tabular}{rrrrrrrrrr}
\toprule
Threads $p$ & $T$ (s) & spread (s) & $S$ & $E$ & $T_P/T$ (\%) & $\bar k$ & $t_{\mathrm{it}}$ (ms) & halo (ms) & sum (ms) \\
\midrule
"""
tail = r"""\bottomrule
\end{tabular}
\caption{Experiment 2: strong scaling over OpenMP threads on the fixed problem $P$ ($1024\times1024$ cells, 100 steps), one rank. $T$ is the median of five repetitions of the loop time and the spread is their largest minus their smallest value; $S=T_{\mathrm{serial}}/T$ and $E=S/p$ are taken against \texttt{heat\_serial}. $T_P/T$ is the share of the Poisson solves, $\bar k$ the mean CG iterations per step and $t_{\mathrm{it}}$, halo and sum the time per CG iteration of the Poisson solve, of its halo routine and of its global sums; with one rank there are no messages, so these two are the time spent inside those routines. The runs are those after the storage was first touched in parallel (8 November 2026, on a GPU node of the same CPUs), which has no one-thread point; the serial reference is that of the first run of the experiment. One socket holds 36 cores, so 64 and 72 threads use both.}
\label{tab:strong-openmp}
\end{table}
"""
out = head + "\n".join(lines) + "\n" + tail
(Path(__file__).resolve().parents[1] / "results/figures/strong_scaling_openmp_table.tex").write_text(out)
print(out)
