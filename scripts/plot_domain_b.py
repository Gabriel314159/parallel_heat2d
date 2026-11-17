#!/usr/bin/env python3
"""Report figure: the 3D benchmark of Model B, as Figure 1 is for 2D: the closed cube with the initial
pigment blob as a cloud of points, the wind against y and its surface profile against z.

Values follow model_b.cpp: blob at (0.3, 0.5, 0.5) L, width 0.08 L, amplitude 1; wind
f_x = (F0/k_y) cos(k_y y) g(z), g(z) = 4 exp[-4(1 - z/Lz)] / (1 - e^-4).
"""
import numpy as np

import fieldviz
import plot_style as ps

L = 1.0
x0, y0, z0, sigma = 0.3, 0.5, 0.5, 0.08

ps.apply()
ps.plt.rcParams.update({"axes.grid": False, "axes.spines.top": True, "axes.spines.right": True,
                        "mathtext.fontset": "cm"})
fig = ps.plt.figure(figsize=(4.6, 4.5))
note = {"fontsize": 10, "color": ps.INK}

# (a) the cube and the blob
ax = fig.add_axes([-0.04, -0.08, 1.04, 1.0], projection="3d", computed_zorder=False)
for z in (0, L):
    ax.plot([0, L, L, 0, 0], [0, 0, L, L, 0], [z] * 5, color=ps.INK, linewidth=0.9)
for x, y in ((0, 0), (L, 0), (L, L), (0, L)):
    ax.plot([x, x], [y, y], [0, L], color=ps.INK, linewidth=0.9)
xx, yy = np.meshgrid([0, L], [0, L])
ax.plot_surface(xx, yy, np.full_like(xx, L), color=ps.BLUE, alpha=0.06, linewidth=0, shade=False)
rng = np.random.default_rng(1)
pts = rng.normal([x0, y0, z0], sigma, size=(6000, 3))
pts = pts[np.all((pts > 0) & (pts < L), axis=1)]
c = np.exp(-np.sum((pts - [x0, y0, z0]) ** 2, axis=1) / (2 * sigma ** 2))
order = np.argsort(c)
norm, cmap = fieldviz.make_norm(0.0, 1.0, 0.5)
ax.scatter(*pts[order].T, c=c[order], cmap=cmap, norm=norm, s=1.6, linewidths=0, depthshade=False)
ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
ax.set_box_aspect((1, 1, 1))
for setter in (ax.set_xticks, ax.set_yticks, ax.set_zticks):
    setter([0, 0.5, 1])
for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis.set_major_formatter("{x:g}")
    axis.pane.set_visible(False)
    axis._axinfo["grid"]["color"] = ps.GRID
ax.tick_params(labelsize=8, pad=0)
ax.set_xlabel("$x$", labelpad=-4); ax.set_ylabel("$y$", labelpad=-4); ax.set_zlabel("$z$", labelpad=-4)
ax.view_init(elev=22, azim=-60)
fig.text(0.5, 0.95, r"$\mathbf{u}\cdot\mathbf{n}=0,\ \ \partial c/\partial n=0\ \ \mathrm{on}\ \partial\Omega$",
         ha="center", **note)

# the wind profiles drawn on the walls: cos(k_y y) on the floor (xy plane), g(z) on the back wall x = 0 (yz plane)
s = np.linspace(0, L, 200)
kw = {"color": ps.ORANGE, "linewidth": 1.4, "linestyle": ":"}
ref = {"color": ps.MUTED, "linewidth": 0.6, "linestyle": "-", "alpha": 0.6}
ax.plot(np.full_like(s, 0.5), s, 0 * s, **ref)
ax.plot(0.5 + 0.3 * np.cos(2 * np.pi * s / L), s, 0 * s, **kw)
g = 4 * np.exp(-4 * (1 - s / L)) / (1 - np.exp(-4))
ax.plot(0 * s, np.full_like(s, 0.2 * 1), s, **ref)
ax.plot(0 * s, 0.2 * g, s, **kw)
ax.text(0.68, 0.3, 0, r"$\cos(k_y y)$", color=ps.ORANGE, fontsize=9)
ax.text(0, 0.3, 0.4, r"$g(z)$", color=ps.ORANGE, fontsize=9)
ps.save(fig, "domain_b")
