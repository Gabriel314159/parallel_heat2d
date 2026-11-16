"""Experiment 5 demo: where the loop time goes in the weak-scaled runs (share of loop, share and ms per iteration of the Poisson solve)."""
import statistics as st
import sys

import plot_style as ps

ps.apply()
runs = ps.load("weak_scaling_rerun.csv")
ranks = sorted(int(k[1:].split("_")[0]) for k in runs if k.startswith("p"))
rows = {p: runs[f"p{p}_t1_2d_blocking"] for p in ranks}


def med(p, key):
    return st.median(float(r[key]) for r in rows[p])


T = {p: med(p, "elapsed_s") for p in ranks}
P = {p: med(p, "poisson_s") for p in ranks}
halo = {p: med(p, "poisson_halo_s") for p in ranks}
MIN = "min" in sys.argv[1:]  # global sums as the minimum over ranks (no waiting for the slowest rank)
tsum = {p: med(p, "poisson_sum_min_s" if MIN else "poisson_sum_s") for p in ranks}
n = {p: float(rows[p][0]["steps"]) * med(p, "poisson_iterations_mean") for p in ranks}
x = range(len(ranks))
labels = [str(p) for p in ranks]
GREY = "#c9c8c2"
comp = {p: P[p] - halo[p] - tsum[p] for p in ranks}


def stack(ax, parts, colours, names, fmt=None):
    bottom = [0.0] * len(ranks)
    for vals, c, name in zip(parts, colours, names):
        ax.bar(x, vals, bottom=bottom, color=c, edgecolor="white", linewidth=1.5, width=0.72, label=name, zorder=3)
        if fmt:
            for i, (v, b) in enumerate(zip(vals, bottom)):
                if v >= 0.07 * (b + v + 1e-12) and v > fmt[1]:
                    ax.text(i, b + v / 2, fmt[0](v), ha="center", va="center", fontsize=8, color="white" if c != ps.YELLOW else ps.INK)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_xticks(list(x), labels)
    ax.grid(axis="x", visible=False)


fig, axes = ps.plt.subplots(2, 2, figsize=(11.5, 6.4), layout="constrained")
ax = axes[0, 0]
stack(ax, [[100 * P[p] / T[p] for p in ranks], [100 * (T[p] - P[p]) / T[p] for p in ranks]], [ps.BLUE, GREY],
      ["Poisson solve", "everything else"], (lambda v: f"{v:.1f}", 50))
ax.set_ylim(0, 100); ax.set_ylabel("share of the loop time (%)"); ax.set_title("Loop time"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
ax = axes[0, 1]
stack(ax, [[P[p] for p in ranks], [T[p] - P[p] for p in ranks]], [ps.BLUE, GREY], ["Poisson solve", "everything else"])
ax.set_ylabel("loop time $T$ (s)"); ax.set_title("Loop time, absolute"); ax.legend(loc="upper left")
parts = [[comp[p] / P[p] * 100 for p in ranks], [halo[p] / P[p] * 100 for p in ranks], [tsum[p] / P[p] * 100 for p in ranks]]
names = ["compute (rest)", "halo exchange", "global sums (min over ranks)" if MIN else "global sums"]
cols = [ps.PURPLE, ps.AQUA, ps.YELLOW]
ax = axes[1, 0]
stack(ax, parts, cols, names, (lambda v: f"{v:.0f}", 0))
ax.set_ylim(0, 100); ax.set_ylabel("share of the Poisson solve (%)"); ax.set_title("Poisson solve"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
ax = axes[1, 1]
stack(ax, [[1e3 * v[p] / n[p] for p in ranks] for v in (comp, halo, tsum)], cols, names)
ax.set_ylabel("time per CG iteration $t_{it}$ (ms)"); ax.set_title("Poisson solve, absolute"); ax.legend(loc="upper left")
for ax in axes[1]:
    ax.set_xlabel("MPI ranks (cores)")
for ax in axes[0]:
    ax.set_xlabel("MPI ranks (cores)")
fig.suptitle("Weak scaling with MPI: breakdown of the loop", fontsize=13)
ps.save(fig, "weak_scaling_mpi_breakdown" + ("_summin" if MIN else ""))
