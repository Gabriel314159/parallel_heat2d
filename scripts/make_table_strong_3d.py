"""LaTeX table of the Model B strong scaling on 1024^3 cells, from results/raw/model_b_3d/strong_p*/raw.csv.
The 1-rank time is extrapolated from its log, as in plot_strong_3d.py."""
import csv
import re
import statistics as st

import plot_style as ps

DIR = ps.RAW / "model_b_3d"
rows = {}
for path in sorted(DIR.glob("strong_p*/raw.csv")):
    runs = [r for r in csv.DictReader(path.open(newline="")) if r["phase"] == "measured"]
    rows[int(runs[0]["ranks"])] = runs
k = int(float(next(iter(rows.values()))[0]["cg_iterations"]))
progress = [(int(a), float(b)) for a, b in
            re.findall(r"iteration=(\d+) solve_elapsed_s=([0-9.]+)", (DIR / "strong_p1/1x1_measured_01.log").read_text())]
(k_end, t_end), (k_prev, t_prev) = progress[-1], progress[-3]
t1 = t_end + (k - k_end) * (t_end - t_prev) / (k_end - k_prev)


def med(runs, key):
    return st.median(float(r[key]) for r in runs)


lines = [r"1 & $1\times1\times1$ & %.0f$^\dagger$ & 1.0 & 1.00 & -- & %.0f & -- & -- \\" % (t1, t1 / k * 1e3)]
for p, runs in sorted(rows.items()):
    T, TP = med(runs, "loop_s"), med(runs, "poisson_s")
    grid = r"$%s\times%s\times%s$" % (runs[0]["px"], runs[0]["py"], runs[0]["pz"])
    lines.append(r"%d & %s & %.0f & %.1f & %.2f & %.1f & %.0f & %.1f & %.1f \\" % (
        p, grid, T, t1 / T, t1 / T / p, 100 * TP / T, TP / k * 1e3,
        med(runs, "halo_s") / k * 1e3, med(runs, "reduction_s") / k * 1e3))

tex = r"""\begin{table}[htbp]
\centering
\small
\begin{tabular}{rcrrrrrrr}
\toprule
Ranks $p$ & Process grid & $T$ (s) & $S$ & $E$ & $T_P/T$ (\%%) & $t_{\mathrm{it}}$ (ms) & halo (ms) & sum (ms) \\
\midrule
%s
\bottomrule
\end{tabular}
\caption{Strong scaling of the three-dimensional model over MPI ranks on the fixed problem $P_3$ ($1024^3$ cells, one time step, %d CG iterations), one thread per rank, on one node. $T$ is the loop time (the median of three repetitions at 64 ranks, a single run otherwise), and $S=T(1)/T$ and $E=S/p$ are taken against the run on one rank. $^\dagger$The one-rank run reached its time limit at CG iteration %d; its time is the elapsed time at that point plus the remaining %d iterations at the rate of the last 200. $T_P/T$ is the share of the pressure solve, and $t_{\mathrm{it}}$, halo and sum the time per CG iteration of the solve, of its halo exchanges and of its global sums.}
\label{tab:strong-3d}
\end{table}
""" % ("\n".join(lines), k, k_end, k - k_end)
out = ps.FIGURES / "strong_scaling_3d_table.tex"
out.write_text(tex)
print(tex)
