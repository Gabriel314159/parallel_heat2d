#!/usr/bin/env python3
"""Report figure: the 3D pigment of Model B as clouds of points, at release and at later times after it.

Reads the final field of model_b --output (CSV with x, y, z, c, all increasing). Points are drawn with a
density proportional to c and coloured with the scale of flow_evolution.

Example:
    python3 scripts/plot_domain_3d.py outputs/modelb3d/b1.csv:0 outputs/modelb3d/b1800.csv:2.5 outputs/modelb3d/b7200.csv:10
"""
import argparse

import numpy as np

import fieldviz
import plot_style as ps

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fields", nargs="+", help="FILE.csv:TAU, one per panel")
parser.add_argument("--lx", type=float, default=1.0)
parser.add_argument("--ly", type=float, default=1.0)
parser.add_argument("--lz", type=float, default=0.5)
parser.add_argument("--stretch", type=float, default=1.0, help="vertical exaggeration (default 1)")
parser.add_argument("--points", type=int, default=9000, help="points per panel")
parser.add_argument("--name", default="domain_3d")
args = parser.parse_args()
lx, ly, lz, s = args.lx, args.ly, args.lz, args.stretch



def load(path):
    """Concentration as an array (nz, ny, nx), indices increasing with z, y, x."""
    d = np.genfromtxt(path, delimiter=",", names=True)
    d = d[np.lexsort((d["x"], d["y"], d["z"]))]
    n = [len(np.unique(d[k])) for k in "zyx"]
    return d["c"].reshape(n)


panels = [(load(f.rsplit(":", 1)[0]), float(f.rsplit(":", 1)[1])) for f in args.fields]
peak = max(float(c.max()) for c, _ in panels)
norm, cmap = fieldviz.make_norm(0.0, peak, 0.5)
rng = np.random.default_rng(1)


def cloud(c):
    """Points spread over the cells with probability proportional to c, jittered inside each cell."""
    nz, ny, nx = c.shape
    flat = np.clip(c, 0, None).ravel()
    cells = rng.choice(flat.size, size=args.points, p=flat / flat.sum())
    k, i, j = np.unravel_index(cells, c.shape)
    x = (j + rng.uniform(size=cells.size)) * lx / nx
    y = (i + rng.uniform(size=cells.size)) * ly / ny
    z = (k + rng.uniform(size=cells.size)) * lz / nz
    order = np.argsort(flat[cells])
    return x[order], y[order], z[order], flat[cells][order]


ps.apply()
ps.plt.rcParams.update({"mathtext.fontset": "cm"})
n = len(panels)
fig = ps.plt.figure(figsize=(4.0 * n + 0.8, 3.0))
width = 0.88 / n
for p, (c, tau) in enumerate(panels):
    ax = fig.add_axes([p * width - 0.01, -0.14, width, 1.24], projection="3d", computed_zorder=False)
    for z in (0, lz):
        ax.plot([0, lx, lx, 0, 0], [0, 0, ly, ly, 0], [z * s] * 5, color=ps.INK, linewidth=0.8)
    for x, y in ((0, 0), (lx, 0), (lx, ly), (0, ly)):
        ax.plot([x, x], [y, y], [0, lz * s], color=ps.INK, linewidth=0.8)
    xx, yy = np.meshgrid([0, lx], [0, ly])
    ax.plot_surface(xx, yy, np.full_like(xx, lz * s), color=ps.BLUE, alpha=0.05, linewidth=0, shade=False)
    x, y, z, v = cloud(c)
    points = ax.scatter(x, y, z * s, c=v, cmap=cmap, norm=norm, s=1.6, linewidths=0, depthshade=False)
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    ax.set_zlim(0, lz * s)
    ax.set_box_aspect((lx, ly, lz * s))
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_zticks([0, lz * s])
    ax.set_zticklabels(["0", f"{lz:g}"] if p == n - 1 else ["", ""])
    ax.xaxis.set_major_formatter("{x:g}")
    ax.yaxis.set_major_formatter("{x:g}")
    ax.tick_params(labelsize=8, pad=0)
    ax.set_xlabel("$x$", labelpad=-4)
    ax.set_ylabel("$y$", labelpad=-4)
    ax.set_zlabel("$z$" if p == n - 1 else "", labelpad=-6)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_visible(False)
        axis._axinfo["grid"]["color"] = ps.GRID
    ax.view_init(elev=22, azim=-60)
    label = "0" if tau == 0 else f"{tau:g}"
    fig.text(p * width + width / 2, 0.9, f"({'abcd'[p]}) $\\tau = {label}$", ha="center", fontsize=11)
bar = fig.colorbar(points, cax=fig.add_axes([0.925, 0.2, 0.014, 0.6]),
                   ticks=[0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0], format="%g")
bar.set_label("pigment concentration $c$")
bar.outline.set_edgecolor(ps.MUTED)
ps.save(fig, args.name)
