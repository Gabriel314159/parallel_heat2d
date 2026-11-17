#!/usr/bin/env python3
"""Report figure: the 3D counterpart of flow_evolution, from Model B runs. (a) The flow at release, as
grey arrows of the horizontal velocity at a few depths, with the blob; (b)-(d) the pigment at three times
as clouds of points whose density follows c.

Reads the final fields of model_b --output (CSV with x, y, z, c, u_e, v_n, w_t, all increasing). The
first file also gives the flow at release.

Example:
    python3 scripts/plot_domain_3d.py outputs/modelb3d/b1.csv:0 outputs/modelb3d/b1800.csv:2.5 \\
        outputs/modelb3d/b7200.csv:10
"""
import argparse

import numpy as np

import fieldviz
import plot_style as ps

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("fields", nargs=3, help="FILE.csv:TAU, one per pigment panel")
parser.add_argument("--lx", type=float, default=1.0)
parser.add_argument("--ly", type=float, default=1.0)
parser.add_argument("--lz", type=float, default=0.5)
parser.add_argument("--points", type=int, default=9000, help="points per panel")
parser.add_argument("--name", default="domain_3d")
args = parser.parse_args()
lx, ly, lz = args.lx, args.ly, args.lz


def load(path):
    """Columns as arrays (nz, ny, nx), indices increasing with z, y, x."""
    d = np.genfromtxt(path, delimiter=",", names=True)
    d = d[np.lexsort((d["x"], d["y"], d["z"]))]
    n = [len(np.unique(d[k])) for k in "zyx"]
    return {k: d[k].reshape(n) for k in d.dtype.names}


shots = [(load(f.rsplit(":", 1)[0]), float(f.rsplit(":", 1)[1])) for f in args.fields]
norm, cmap = fieldviz.make_norm(0.0, float(shots[0][0]["c"].max()), 0.5)
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


def box(ax):
    for z in (0, lz):
        ax.plot([0, lx, lx, 0, 0], [0, 0, ly, ly, 0], [z] * 5, color=ps.INK, linewidth=0.8)
    for x, y in ((0, 0), (lx, 0), (lx, ly), (0, ly)):
        ax.plot([x, x], [y, y], [0, lz], color=ps.INK, linewidth=0.8)
    xx, yy = np.meshgrid([0, lx], [0, ly])
    ax.plot_surface(xx, yy, np.full_like(xx, lz), color=ps.BLUE, alpha=0.05, linewidth=0, shade=False)
    ax.set_xlim(0, lx); ax.set_ylim(0, ly); ax.set_zlim(0, lz)
    ax.set_box_aspect((lx, ly, lz), zoom=0.82)
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.set_zticks([0, lz])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_formatter("{x:g}")
        axis.pane.set_visible(False)
        axis._axinfo["grid"]["color"] = ps.GRID
    ax.tick_params(labelsize=8, pad=0)
    ax.set_xlabel("$x$", labelpad=-4); ax.set_ylabel("$y$", labelpad=-4); ax.set_zlabel("$z$", labelpad=-6)
    ax.view_init(elev=22, azim=-60)


ps.apply()
ps.plt.rcParams.update({"mathtext.fontset": "cm"})
fig = ps.plt.figure(figsize=(11.5, 3.2))
width = 0.222
axes = [fig.add_axes([p * width, -0.04, width, 1.0], projection="3d", computed_zorder=False)
        for p in range(4)]

# (a) the flow at release: streamlines of the depth-averaged horizontal flow, drawn on the surface, and the blob
ax = axes[0]
box(ax)
f = shots[0][0]
nz, ny, nx = f["c"].shape
ubar, vbar = f["u_e"].mean(axis=0), f["v_n"].mean(axis=0)  # (ny, nx): east and north faces
psi = np.vstack([np.zeros(nx), np.cumsum(ubar, axis=0) * ly / ny])  # psi(y_face, x_east), zero on y = 0
xs, ys = (np.arange(nx) + 1) * lx / nx, np.arange(ny + 1) * ly / ny
span = np.abs(psi).max()
scratch = ps.plt.figure()
lines = scratch.add_subplot().contour(xs, ys, psi, levels=np.linspace(-span, span, 10)[1:-1])
ps.plt.close(scratch)
for level, segments in zip(lines.levels, lines.allsegs):
    for path in segments:
        ax.plot(path[:, 0], path[:, 1], np.full(len(path), lz), color=ps.MUTED, linewidth=0.5, alpha=0.7)
    if abs(abs(level) - 0.5 * span) > 0.12 * span or not segments:  # arrows on one contour per gyre
        continue
    path = max(segments, key=len)
    for m in (len(path) // 4, 3 * len(path) // 4):
        (x1, y1), (x2, y2) = path[m], path[m + 1]
        i, j = min(int(x1 / lx * nx), nx - 1), min(int(y1 / ly * ny), ny - 1)
        if (x2 - x1) * ubar[j, i] + (y2 - y1) * vbar[j, i] < 0:
            (x1, y1), (x2, y2) = (x2, y2), (x1, y1)
        ax.quiver(x1, y1, lz, x2 - x1, y2 - y1, 0, length=0.07, normalize=True, color=ps.MUTED,
                  linewidth=0.8, arrow_length_ratio=1.0)
x, y, z, _ = cloud(f["c"])
ax.scatter(x, y, z, color=ps.BLUE, alpha=0.08, s=1.2, linewidths=0, depthshade=False)
ax.set_title("(a) Flow at release", fontsize=11, y=0.98)

for ax, (field, tau), letter in zip(axes[1:], shots, "bcd"):
    box(ax)
    x, y, z, v = cloud(field["c"])
    points = ax.scatter(x, y, z, c=v, cmap=cmap, norm=norm, s=1.6, linewidths=0, depthshade=False)
    ax.set_title(f"({letter}) $\\tau = {tau:g}$", fontsize=11, y=0.98)

bar = fig.colorbar(points, cax=fig.add_axes([0.905, 0.22, 0.011, 0.56]),
                   ticks=[0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0], format="%g")
bar.set_label("pigment concentration $c$")
bar.outline.set_edgecolor(ps.MUTED)
ps.save(fig, args.name)
for field, tau in shots:
    print(f"tau = {tau:g}: max c = {field['c'].max():.3f}")
