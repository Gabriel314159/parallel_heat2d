"""LaTeX table of experiment 5 (weak scaling over MPI ranks) from results/raw."""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
groups = defaultdict(list)
for r in csv.DictReader((ROOT / "results/raw/weak_scaling_rerun.csv").open(newline="")):
    if r["status"] == "ok" and r["phase"] == "measured":
        groups[r["run_id"]].append(r)


def med(rows, key):
    return st.median(float(r[key]) for r in rows)


def stats(rows):
    t = [float(r["elapsed_s"]) for r in rows]
    n = float(rows[0]["steps"]) * med(rows, "poisson_iterations_mean")
    return st.median(t), max(t) - min(t), med(rows, "poisson_iterations_mean"), {
        c: 1e3 * med(rows, c) / n for c in ("poisson_s", "poisson_halo_s", "poisson_sum_s", "poisson_sum_min_s")}


T1, _, _, one = stats(groups["p1_t1_2d_blocking"])
lines = []
t = [float(r["elapsed_s"]) for r in groups["serial_reference"]]
lines.append(f"serial & $256^2$ & {st.median(t):.2f} & " + " & ".join(["--"] * 6) + r" \\")
lines.append(r"\midrule")
for p in (1, 4, 9, 16, 36, 64, 144, 256):
    rows = groups[f"p{p}_t1_2d_blocking"]
    T, spread, k, c = stats(rows)
    nx = int(rows[0]["nx"])
    lines.append(f"{p} & ${nx}^2$ & {T:.2f} & {T1 / T:.3f} & {k:.0f} & {c['poisson_s']:.3f} & "
                 f"{one['poisson_s'] / c['poisson_s']:.2f} & {c['poisson_sum_s']:.3f} & {c['poisson_sum_min_s']:.3f} \\\\")
head = r"""\begin{table}[!t]
\centering
\small
\begin{tabular}{rrrrrrrrr}
\toprule
Ranks $p$ & Grid & $T$ (s) & $E_w$ & $\bar k$ & $t_{\mathrm{it}}$ (ms) & $E_{w,\mathrm{it}}$ & sum (ms) & min sum (ms) \\
\midrule
"""
tail = r"""\bottomrule
\end{tabular}
\caption{Experiment 5: weak scaling over MPI ranks, medians of five repetitions. $E_w=T(1)/T(p)$ and $E_{w,\mathrm{it}}=t_{\mathrm{it}}(1)/t_{\mathrm{it}}(p)$ are taken against the one-rank run. Sum is the time per CG iteration inside the global sums for the slowest rank and min sum for the fastest rank. Up to 64 ranks the runs are on one node; 144 ranks use three nodes and 256 ranks four.}
\label{tab:weak-mpi}
\end{table}
"""
out = head + "\n".join(lines) + "\n" + tail
(ROOT / "results/figures/weak_scaling_mpi_table.tex").write_text(out)
print(out)
