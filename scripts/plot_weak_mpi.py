"""Experiment 5, weak scaling over MPI ranks: runtime, time per CG iteration, and the two efficiencies."""
import statistics as st

import plot_style as ps

ps.apply()
runs = ps.load("weak_scaling_rerun.csv")
ranks = sorted(int(k[1:].split("_")[0]) for k in runs if k.startswith("p"))
rows = {p: runs[f"p{p}_t1_2d_blocking"] for p in ranks}
T = {p: ps.stats(rows[p]) for p in ranks}
k = {p: st.median(float(r["poisson_iterations_mean"]) for r in rows[p]) for p in ranks}
tit = {p: ps.per_iteration(rows[p]) * 1e3 for p in ranks}
TICKS = [1, 4, 16, 64, 256]

fig, axes = ps.plt.subplots(1, 3, figsize=(11.5, 3.6), layout="constrained")
ax = axes[0]
ps.curve(ax, ranks, [T[p][0] for p in ranks], "serial", "measured", lo=[T[p][1] for p in ranks], hi=[T[p][2] for p in ranks])
ps.curve(ax, ranks, [T[1][0]] * len(ranks), "ideal", "ideal: $T(1)$")
ps.curve(ax, ranks, [T[1][0] * k[p] / k[1] for p in ranks], "ideal_k", "ideal $t_{it}$: $T(1)\\,\\bar{k}(p)/\\bar{k}(1)$")
ax.set_yscale("log"); ax.set_ylabel("runtime $T$ (s)"); ax.set_title("Runtime"); ax.legend(loc="upper left")
ax = axes[1]
ps.curve(ax, ranks, [tit[p] for p in ranks], "serial", "measured")
ps.curve(ax, ranks, [tit[1]] * len(ranks), "ideal", "ideal: $t_{it}(1)$")
ax.set_ylim(0, 0.5); ax.set_ylabel("time per CG iteration $t_{it}$ (ms)"); ax.set_title("Time per iteration")
ax.legend(loc="upper left")
ax = axes[2]
ps.curve(ax, ranks, [T[1][0] / T[p][0] for p in ranks], "serial", "$E_w=T(1)/T(p)$")
ps.curve(ax, ranks, [tit[1] / tit[p] for p in ranks], "serial_omp", "$E_{w,it}=t_{it}(1)/t_{it}(p)$")
ps.curve(ax, ranks, [1] * len(ranks), "ideal", "ideal: $E=1$")
ax.set_ylim(0, 1.7); ax.set_ylabel("weak efficiency"); ax.set_title("Efficiency"); ax.legend(loc="upper right")
for ax in axes:
    ps.ranks_axis(ax, TICKS)
fig.suptitle("Weak scaling with MPI", fontsize=13)
ps.save(fig, "weak_scaling_mpi")
