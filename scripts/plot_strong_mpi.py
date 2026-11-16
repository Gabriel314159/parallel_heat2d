"""Experiment 1, strong scaling over MPI ranks: runtime, speedup, efficiency, and the time per CG iteration."""
import plot_style as ps

ps.apply()
runs = ps.load("strong_scaling_mpi_1node.csv", "strong_scaling_mpi_4nodes.csv")
ranks = sorted(int(k[1:].split("_")[0]) for k in runs if k.startswith("p"))
rows = {p: runs[f"p{p}_t1_2d_blocking"] for p in ranks}
ref = {"serial": ps.stats(runs["serial_reference"])[0], "serial_omp": ps.stats(runs["serial_omp_reference"])[0]}
T = {p: ps.stats(rows[p]) for p in ranks}

fig, axes = ps.plt.subplots(1, 3, figsize=(11.5, 3.6), layout="constrained")
ax = axes[0]
ps.curve(ax, ranks, [T[p][0] for p in ranks], "serial", "measured", lo=[T[p][1] for p in ranks], hi=[T[p][2] for p in ranks])
ps.curve(ax, ranks, [ref["serial"] / p for p in ranks], "ideal", "ideal: $T_{serial}/p$")
ax.set_yscale("log"); ax.set_ylabel("runtime $T$ (s)"); ax.set_title("Runtime")
ax.legend()
ax = axes[1]
ps.curve(ax, ranks, [ref["serial"] / T[p][0] for p in ranks], "serial", "measured")
ps.curve(ax, ranks, ranks, "ideal", "ideal: $S=p$")
ax.set_yscale("log", base=2); ax.set_ylabel("speedup $S(p)$"); ax.set_title("Speedup")
ax.yaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v:g}"))
ax.legend()
ax = axes[2]
ps.curve(ax, ranks, [ref["serial"] / T[p][0] / p for p in ranks], "serial", "measured")
ps.curve(ax, ranks, [1] * len(ranks), "ideal", "ideal: $E=1$")
ax.legend(loc="lower left")
ax.set_ylim(0, 1.25); ax.set_ylabel("efficiency $E(p)=S(p)/p$"); ax.set_title("Efficiency")
for ax in axes:
    ps.ranks_axis(ax, ranks[::2] + ([ranks[-1]] if len(ranks) % 2 == 0 else []))
fig.suptitle("Strong scaling with MPI", fontsize=13)
ps.save(fig, "strong_scaling_mpi")

fig, ax = ps.plt.subplots(figsize=(5.5, 3.8))
t = {p: ps.per_iteration(rows[p]) * 1e3 for p in ranks}
ps.curve(ax, ranks, [t[p] for p in ranks], "serial", "measured")
ps.curve(ax, ranks, [t[1] / p for p in ranks], "ideal", "ideal: $t_{it}(1)/p$")
ax.set_yscale("log"); ax.set_ylabel("time per CG iteration (ms)")
ax.set_title("Strong scaling with MPI")
ps.ranks_axis(ax, ranks[::2] + ([ranks[-1]] if len(ranks) % 2 == 0 else []))
ax.legend()
ps.save(fig, "strong_scaling_mpi_iteration")
