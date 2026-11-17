"""Cost of one CG iteration per cell on one core, against the CG working set of a rank, with the
cache sizes of the Xeon Platinum 8452Y marked. Points: serial runs of growing grids (2D), the compute
part of MPI strong scaling on 1024^2 (2D, halo and sums removed) and on 1024^3 (3D, Model B).

The working set is the four CG arrays a rank sweeps every iteration (psi or phi, r, p, Ap), 8 B per cell.
"""
import csv
import statistics as st

import plot_style as ps

ps.apply()
RAW = ps.RAW
ARRAYS, BYTES = 4, 8


def rows(name):
    return [r for r in csv.DictReader((RAW / name).open(newline="")) if r["phase"] == "measured" and r["status"] == "ok"]


# 2D serial: ns per unknown and iteration, unknowns (n-1)^2
serial = {}
for name in ("serial_sizes.csv", "serial_sizes_2048.csv", "serial_sizes_4096.csv"):
    for r in rows(name):
        if r["backend"] == "serial":
            serial.setdefault(int(r["nx"]), []).append(r)
s2d, iters = [], []
for n, rs in sorted(serial.items()):
    unknowns = (n - 1) ** 2
    t_it = st.median(float(r["poisson_s"]) / (int(r["steps"]) * float(r["poisson_iterations_mean"])) for r in rs)
    s2d.append((n, unknowns * ARRAYS * BYTES, t_it / unknowns * 1e9))
    iters.append((n, st.median(float(r["poisson_iterations_mean"]) for r in rs)))

# 2D strong scaling on one node: compute part of an iteration
strong = {}
for r in rows("strong_scaling_mpi_1node.csv"):
    if r["backend"] != "serial":
        strong.setdefault(int(r["mpi_ranks"]), []).append(r)
m2d = []
for p, rs in sorted(strong.items()):
    unknowns = 1023 ** 2
    t = st.median((float(r["poisson_s"]) - float(r["poisson_halo_s"]) - float(r["poisson_sum_s"]))
                  / (int(r["steps"]) * float(r["poisson_iterations_mean"])) for r in rs)
    m2d.append((p, unknowns / p * ARRAYS * BYTES, t * p / unknowns * 1e9))

# 3D strong scaling (Model B), 1024^3, compute part; the one-rank run is the extrapolated point
m3d = []
for path in sorted((RAW / "model_b_3d").glob("strong_p*/raw.csv")):
    rs = [r for r in csv.DictReader(path.open(newline="")) if r["phase"] == "measured"]
    p, k, cells = int(rs[0]["ranks"]), float(rs[0]["cg_iterations"]), 1024 ** 3
    t = st.median((float(r["poisson_s"]) - float(r["halo_s"]) - float(r["reduction_s"])) / k for r in rs)
    m3d.append((p, cells / p * ARRAYS * BYTES, t * p / cells * 1e9))
m3d.sort()

MB = 2 ** 20
fig, (left, ax) = ps.plt.subplots(1, 2, figsize=(10.5, 4.0), layout="constrained", gridspec_kw={"width_ratios": [1, 1.6]})

# (a) CG iterations per step against n, 2D serial, with a slope-one reference through 1024
ns, ks = [n for n, _ in iters], [k for _, k in iters]
k0 = dict(iters)[1024]
left.plot(ns, [k0 * n / 1024 for n in ns], "--", color=ps.MUTED, linewidth=1, label=r"$\propto n$")
left.plot(ns, ks, "o-", color=ps.MUTED, label="2D serial")
for n, k in iters:
    left.annotate(f"{k:.0f}", (n, k), textcoords="offset points", xytext=(-4, 7), ha="right", fontsize=8, color=ps.MUTED)
left.set_xscale("log", base=2); left.set_yscale("log")
left.set_xticks(ns); left.xaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"${int(v)}^2$"))
left.yaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v:g}"))
left.set_xlabel(r"grid $n\times n$"); left.set_ylabel("CG iterations per step")
left.set_title("(a) Iterations", fontsize=11, loc="left")
left.legend(loc="upper left", fontsize=8)
ax.set_title("(b) Cost per cell and iteration", fontsize=11, loc="left")
L2, L3 = 2 * MB, 67.5 * MB  # per core; per socket, shared by up to 36 cores
ax.axvspan(2.5e5, L2, color=ps.BLUE, alpha=0.07, linewidth=0)
ax.axvspan(L2, L3, color=ps.BLUE, alpha=0.03, linewidth=0)
for x, label in ((L2, "L2\n2 MB per core"), (L3, "L3\n67.5 MB per socket")):
    ax.axvline(x, color=ps.MUTED, linewidth=0.8, linestyle="--")
    ax.text(x / 1.1, 19, label, color=ps.MUTED, fontsize=8, va="top", ha="right")

ax.plot([b for _, b, _ in s2d], [c for *_, c in s2d], "o-", color=ps.MUTED, label=r"2D serial, $n^2$ grids")
for n, b, c in s2d:
    ax.annotate(f"${n}^2$", (b, c), textcoords="offset points", xytext=(0, -12), ha="center", fontsize=8, color=ps.MUTED)
ax.plot([b for _, b, _ in m2d], [c for *_, c in m2d], "s-", color=ps.BLUE, label=r"2D MPI, $1024^2$, compute only")
for p, b, c in m2d:
    ax.annotate(f"{p}", (b, c), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8, color=ps.BLUE)
ax.plot([b for _, b, _ in m3d], [c for *_, c in m3d], "^-", color=ps.ORANGE, label=r"3D MPI, $1024^3$, compute only")
for p, b, c in m3d:
    ax.annotate(f"{p}", (b, c), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8, color=ps.ORANGE)

ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlim(2.5e5, 6e10)
ax.set_ylim(1.2, 22)
ax.xaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v / MB:g} MB" if v < 1024 * MB else f"{v / 1024 / MB:g} GB"))
ax.yaxis.set_major_formatter(ps.FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_xlabel("CG working set per rank (4 arrays)")
ax.set_ylabel("ns per cell and CG iteration, one core")
ax.legend(loc="center right", fontsize=8)
ps.save(fig, "cell_cost")
for name, pts in (("2D serial", s2d), ("2D MPI", m2d), ("3D MPI", m3d)):
    print(name, [(a, round(b / MB, 2), round(c, 2)) for a, b, c in pts])
