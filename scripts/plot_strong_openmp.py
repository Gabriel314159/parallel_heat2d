"""Experiment 2, strong scaling over OpenMP threads: runtime, speedup, efficiency, and the time per CG iteration."""
import plot_style as ps

ps.apply()
ref = ps.load("strong_scaling_openmp.csv")["serial_reference"]
runs = ps.load("strong_scaling_openmp_firsttouch.csv")
t_ser = ps.stats(ref)[0]
SETS = (("measured", "serial", runs),)
threads = {name: sorted(int(k.split("_t")[1].split("_")[0]) for k in runs if k.startswith("p1_t")) for name, _, runs in SETS}
TICKS = [1, 4, 16, 36, 72]


def rows(runs, t):
    return runs[f"p1_t{t}_2d_blocking"]


fig, axes = ps.plt.subplots(1, 3, figsize=(11.5, 3.6), layout="constrained")
for name, role, runs in SETS:
    x = threads[name]
    T = [ps.stats(rows(runs, t)) for t in x]
    ps.curve(axes[0], x, [v[0] for v in T], role, name, lo=[v[1] for v in T], hi=[v[2] for v in T])
    ps.curve(axes[1], x, [t_ser / v[0] for v in T], role, name)
    ps.curve(axes[2], x, [t_ser / v[0] / t for v, t in zip(T, x)], role, name)
x = threads["measured"]
ps.curve(axes[0], x, [t_ser / t for t in x], "ideal", "ideal: $T_{serial}/p$")
ps.curve(axes[1], x, x, "ideal", "ideal: $S=p$")
ps.curve(axes[2], x, [1] * len(x), "ideal", "ideal: $E=1$")
axes[0].set_yscale("log"); axes[0].set_ylabel("runtime $T$ (s)"); axes[0].set_title("Runtime")
axes[1].set_yscale("log", base=2); axes[1].set_ylabel("speedup $S(p)$"); axes[1].set_title("Speedup")
axes[1].yaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v:g}"))
axes[2].set_ylim(0, 1.25); axes[2].set_ylabel("efficiency $E(p)=S(p)/p$"); axes[2].set_title("Efficiency")
for ax in axes:
    ps.ranks_axis(ax, TICKS, "OpenMP threads (cores)")
axes[0].legend(); axes[1].legend(loc="upper left"); axes[2].legend(loc="lower left")
fig.suptitle("Strong scaling with OpenMP", fontsize=13)
ps.save(fig, "strong_scaling_openmp")

fig, ax = ps.plt.subplots(figsize=(5.5, 3.8))
for name, role, runs in SETS:
    x = threads[name]
    ps.curve(ax, x, [ps.per_iteration(rows(runs, t)) * 1e3 for t in x], role, name)
t1 = ps.per_iteration(ref) * 1e3
ps.curve(ax, threads["measured"], [t1 / t for t in threads["measured"]], "ideal", "ideal: $t_{it,serial}/p$")
ax.set_yscale("log"); ax.set_ylabel("time per CG iteration (ms)"); ax.set_title("Strong scaling with OpenMP")
ps.ranks_axis(ax, TICKS, "OpenMP threads (cores)")
ax.legend()
ps.save(fig, "strong_scaling_openmp_iteration")
