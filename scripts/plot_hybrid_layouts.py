"""Experiment 3, 256 cores split between ranks and threads: the time per CG iteration and its two parts."""
import plot_style as ps

ps.apply()
runs = ps.load("hybrid_layouts_firsttouch.csv")
layouts = sorted(((int(r[0]["mpi_ranks"]), int(r[0]["omp_threads"])) for r in runs.values()), key=lambda k: k[1])
rows = {t: runs[f"p{p}_t{t}_2d_blocking"] for p, t in layouts}
threads = [t for _, t in layouts]


def per_rep(r, key):
    return 1e3 * float(r[key]) / (float(r["steps"]) * float(r["poisson_iterations_mean"]))


def band(t, f):
    values = [f(r) for r in rows[t]]
    return ps.st.median(values), min(values), max(values)


it = {t: band(t, lambda r: per_rep(r, "poisson_s")) for t in threads}
tsum = {t: band(t, lambda r: per_rep(r, "poisson_sum_s")) for t in threads}
rest = {t: band(t, lambda r: per_rep(r, "poisson_s") - per_rep(r, "poisson_sum_s")) for t in threads}

fig, ax = ps.plt.subplots(figsize=(5.5, 3.8), layout="constrained")
for data, role, label in ((it, "serial", r"$t_{it}$"), (tsum, "sum", r"$t_{sums}$"),
                          (rest, "rest", r"rest: $t_{it}-t_{sums}$")):
    ps.curve(ax, threads, [data[t][0] for t in threads], role, label,
             lo=[data[t][1] for t in threads], hi=[data[t][2] for t in threads])
ax.set_ylim(0, 0.3)  # linear: the two parts add up to t_it
ax.set_ylabel("time per CG iteration (ms)")
ax.set_xscale("log", base=2)
ax.set_xticks(threads)
ax.set_xticklabels([f"{p}×{t}" for p, t in layouts], fontsize=8.5)
ax.xaxis.set_minor_formatter(ps.NullFormatter())
ax.set_xlabel("ranks × threads (256 cores, 4 nodes)")
ax.set_title("Pure MPI vs MPI + OpenMP")
ax.legend(loc="center left")
ps.save(fig, "hybrid_layouts_iteration")
