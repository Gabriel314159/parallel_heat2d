"""Experiment 5, weak scaling over MPI ranks against the serial program on the same grid:
runtime, scaled speedup S_s(p) = T_serial(N_p) / T(p) and scaled efficiency S_s(p) / p."""
import csv
import statistics as st
from collections import defaultdict

import plot_style as ps

ps.apply()
runs = ps.load("weak_scaling_rerun.csv")
ranks = sorted(int(k[1:].split("_")[0]) for k in runs if k.startswith("p"))
rows = {p: runs[f"p{p}_t1_2d_blocking"] for p in ranks}
T = {p: ps.stats(rows[p]) for p in ranks}
tit = {p: ps.per_iteration(rows[p]) for p in ranks}
grid = {p: int(rows[p][0]["nx"]) for p in ranks}

# heat_serial on each weak-scaled grid, by cells a side; the 256^2 one is the rerun's own reference.
serial = defaultdict(list)
for name in ["weak_scaling_rerun.csv"] + sorted(f.name for f in ps.RAW.glob("serial_weak_*.csv")):
    for r in csv.DictReader((ps.RAW / name).open(newline="")):
        if r["status"] == "ok" and r["phase"] == "measured" and r["run_id"] == "serial_reference":
            serial[int(r["nx"])].append(r)
Ts = {p: ps.stats(serial[grid[p]]) for p in ranks}
S = {p: Ts[p][0] / T[p][0] for p in ranks}
TICKS = [1, 4, 16, 64, 256]

fig, axes = ps.plt.subplots(1, 3, figsize=(11.5, 3.6), layout="constrained")
ax = axes[0]
ps.curve(ax, ranks, [T[p][0] for p in ranks], "serial", "measured", lo=[T[p][1] for p in ranks], hi=[T[p][2] for p in ranks])
ps.curve(ax, ranks, [Ts[p][0] for p in ranks], "ideal_k", "serial, same grid: $T_{serial}(N_p)$")
ps.curve(ax, ranks, [Ts[p][0] / p for p in ranks], "ideal", "ideal: $T_{serial}(N_p)/p$")
ax.set_yscale("log"); ax.set_ylabel("runtime $T$ (s)"); ax.set_title("Runtime"); ax.legend(loc="upper left")
ax = axes[1]
ps.curve(ax, ranks, [S[p] for p in ranks], "serial", "measured")
ps.curve(ax, ranks, ranks, "ideal", "ideal: $S_s=p$")
ax.set_yscale("log", base=2); ax.set_ylabel("scaled speedup $S_s(p)$"); ax.set_title("Speedup")
ax.yaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v:g}"))
ax.legend()
ax = axes[2]
ps.curve(ax, ranks, [S[p] / p for p in ranks], "serial", "$E_s=S_s(p)/p$")
ps.curve(ax, ranks, [tit[1] / tit[p] for p in ranks], "serial_omp", "$E_{w,it}=t_{it}(1)/t_{it}(p)$")
ps.curve(ax, ranks, [1] * len(ranks), "ideal", "ideal: $E=1$")
ax.set_ylim(0, 1.7); ax.set_ylabel("efficiency"); ax.set_title("Efficiency"); ax.legend(loc="lower left")
for ax in axes:
    ps.ranks_axis(ax, TICKS)
    ps.mark_nodes(ax, label=ax is axes[1])
fig.suptitle("Weak scaling with MPI", fontsize=13)
ps.save(fig, "weak_scaling_mpi_scaled")

print("p grid T T_serial reps S_s E_s E_w,it E_s/E_w,it k_par k_serial")
for p in ranks:
    ks = st.median(float(r["poisson_iterations_mean"]) for r in serial[grid[p]])
    kp = st.median(float(r["poisson_iterations_mean"]) for r in rows[p])
    e = tit[1] / tit[p]
    print(p, grid[p], f"{T[p][0]:.2f} {Ts[p][0]:.1f} {len(serial[grid[p]])} {S[p]:.1f} {S[p] / p:.2f} {e:.2f} {S[p] / p / e:.2f} {kp:.0f} {ks:.0f}")
