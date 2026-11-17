#!/usr/bin/env python3
"""Report figure: the basin with its walls and the initial pigment blob, with the wind profile
sin(2 pi y / Ly) drawn inside it, in the style of plot_domain_b.py.

Reads the benchmark problem from config/defaults.json; the blob is drawn with the colours of flow_evolution.
"""
import json

import numpy as np

import fieldviz
import plot_style as ps

problem = json.loads((ps.REPO / "config/defaults.json").read_text())["problem"]
lx, ly = problem["lx"], problem["ly"]
x0, y0 = problem["blob_centre"]
sigma, a0 = problem["blob_width"], problem["amplitude"]

ps.apply()
ps.plt.rcParams.update({"axes.grid": False, "axes.spines.top": True, "axes.spines.right": True, "mathtext.fontset": "cm"})
fig, ax = ps.plt.subplots(figsize=(4.2, 4.2), layout="constrained")

# the blob as a cloud of points whose density follows c, as in plot_domain_b.py
rng = np.random.default_rng(1)
pts = rng.normal([x0, y0], sigma, size=(3000, 2))
pts = pts[np.all((pts > 0) & (pts < [lx, ly]), axis=1)]
c = a0 * np.exp(-np.sum((pts - [x0, y0]) ** 2, axis=1) / (2 * sigma ** 2))
order = np.argsort(c)
norm, cmap = fieldviz.make_norm(0.0, a0, 0.5)
ax.scatter(*pts[order].T, c=c[order], cmap=cmap, norm=norm, s=2.0, linewidths=0)
for spine in ax.spines.values():  # the four walls
    spine.set_color(ps.INK)
    spine.set_linewidth(1.6)

# the wind profile F / F0 = sin(2 pi y / Ly), drawn about the line x = lx / 2
s = np.linspace(0, ly, 200)
ax.plot(np.full_like(s, lx / 2), s, color=ps.MUTED, linewidth=0.6, alpha=0.6)
ax.plot(lx / 2 + 0.3 * lx * np.sin(2 * np.pi * s / ly), s, color=ps.BLUE, linewidth=1.4, linestyle=":")
ax.text(0.56 * lx, 0.9 * ly, r"$\sin(2\pi y/L_y)$", color=ps.BLUE, fontsize=14)

ax.set_xlim(0, lx)
ax.set_ylim(0, ly)
ax.set_aspect("equal")
ax.set_xticks([0, lx / 2, lx])
ax.set_yticks([0, ly / 2, ly])
ax.xaxis.set_major_formatter("{x:g}")
ax.yaxis.set_major_formatter("{x:g}")
ax.tick_params(labelsize=12)
ax.set_xlabel("$x$", fontsize=15)
ax.set_ylabel("$y$", fontsize=15)
ps.save(fig, "domain")
