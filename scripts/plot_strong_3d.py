"""Model B (3D flow with pressure), strong scaling over MPI ranks on 1024^3 cells: runtime, speedup,
efficiency and the time per CG iteration. Reads results/raw/model_b_3d/strong_p*/raw.csv, whose columns
are not those of the 2D runs, so ps.load is not used. The 1-rank run hit its time limit before the end
of its only solve: its time is extrapolated from the progress lines of its log (see that folder's README)."""
import csv
import re
import statistics as st

import plot_style as ps

ps.apply()
DIR = ps.RAW / "model_b_3d"

runs = {}
for path in sorted(DIR.glob("strong_p*/raw.csv")):
    rows = [r for r in csv.DictReader(path.open(newline="")) if r["phase"] == "measured"]
    runs[int(rows[0]["ranks"])] = rows
iterations = int(float(next(iter(runs.values()))[0]["cg_iterations"]))  # the same in every layout
T = {p: (st.median(float(r["loop_s"]) for r in rows), min(float(r["loop_s"]) for r in rows),
         max(float(r["loop_s"]) for r in rows)) for p, rows in runs.items()}
t_it = {p: st.median(float(r["poisson_s"]) / float(r["cg_iterations"]) for r in rows) for p, rows in runs.items()}

# The 1-rank run: elapsed time at the last progress line plus the missing iterations at the rate of its last 200.
progress = [(int(m[0]), float(m[1])) for m in
            re.findall(r"iteration=(\d+) solve_elapsed_s=([0-9.]+)", (DIR / "strong_p1/1x1_measured_01.log").read_text())]
(k_end, t_end), (k_prev, t_prev) = progress[-1], progress[-3]
t1 = t_end + (iterations - k_end) * (t_end - t_prev) / (k_end - k_prev)
t1_it = t1 / iterations

ranks = sorted(runs)
all_ranks = [1] + ranks
T_all = {1: t1, **{p: T[p][0] for p in ranks}}
it_all = {1: t1_it, **t_it}
ticks = all_ranks[::2] + ([all_ranks[-1]] if len(all_ranks) % 2 == 0 else [])
OPEN = dict(markerfacecolor="white", markeredgecolor=ps.BLUE, markeredgewidth=1.6, linestyle="none", zorder=4)


def measured(ax, y_all, lo=None, hi=None):
    """Measured ranks as the solid curve, the extrapolated 1-rank point hollow."""
    ps.curve(ax, ranks, [y_all[p] for p in ranks], "serial", "measured", lo=lo, hi=hi)
    ax.plot([1], [y_all[1]], marker="o", markersize=7, label="1 rank, extrapolated", **OPEN)


def node_line(ax):
    if all_ranks[-1] > 64:
        ps.mark_nodes(ax)


fig, axes = ps.plt.subplots(1, 3, figsize=(11.5, 3.6), layout="constrained")
ax = axes[0]
measured(ax, T_all, lo=[T[p][1] for p in ranks], hi=[T[p][2] for p in ranks])
ps.curve(ax, all_ranks, [t1 / p for p in all_ranks], "ideal", "ideal: $T_1/p$")
ax.set_yscale("log"); ax.set_ylabel("runtime $T$ (s)"); ax.set_title("Runtime")
ax.legend()
ax = axes[1]
measured(ax, {p: t1 / T_all[p] for p in all_ranks})
ps.curve(ax, all_ranks, all_ranks, "ideal", "ideal: $S=p$")
ax.set_yscale("log", base=2); ax.set_ylabel("speedup $S(p)$"); ax.set_title("Speedup")
ax.yaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v:g}"))
ax.legend()
ax = axes[2]
measured(ax, {p: t1 / T_all[p] / p for p in all_ranks})
ps.curve(ax, all_ranks, [1] * len(all_ranks), "ideal", "ideal: $E=1$")
ax.legend(loc="lower left")
ax.set_ylim(0, 1.25); ax.set_ylabel("efficiency $E(p)=S(p)/p$"); ax.set_title("Efficiency")
for ax in axes:
    ps.ranks_axis(ax, ticks)
    node_line(ax)
fig.suptitle("Strong scaling with MPI, 3D model", fontsize=13)
ps.save(fig, "strong_scaling_3d")

fig, ax = ps.plt.subplots(figsize=(5.5, 3.8))
measured(ax, {p: it_all[p] * 1e3 for p in all_ranks})
ps.curve(ax, all_ranks, [it_all[1] * 1e3 / p for p in all_ranks], "ideal", "ideal: $t_{it}(1)/p$")
ax.set_yscale("log"); ax.set_ylabel("time per CG iteration (ms)")
ax.set_title("Strong scaling with MPI, 3D model")
ps.ranks_axis(ax, ticks)
node_line(ax)
ax.legend()
ps.save(fig, "strong_scaling_3d_iteration")
