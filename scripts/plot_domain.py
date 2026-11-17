#!/usr/bin/env python3
"""Report figure: the basin with its walls and the initial pigment blob, and the wind forcing against y.

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
fig, (ax, side) = ps.plt.subplots(1, 2, figsize=(5.5, 3.9), layout="constrained", sharey=True,
                                  gridspec_kw={"width_ratios": [3.0, 1.0]})

n = 400
x = (np.arange(n) + 0.5) * lx / n
y = (np.arange(n) + 0.5) * ly / n
blob = a0 * np.exp(-((x[None, :] - x0) ** 2 + (y[::-1, None] - y0) ** 2) / (2 * sigma ** 2))
norm, cmap = fieldviz.make_norm(0.0, a0, 0.5)
ax.imshow(blob, extent=(0, lx, 0, ly), origin="upper", cmap=cmap, norm=norm, interpolation="bilinear")
for spine in ax.spines.values():  # the four walls
    spine.set_color(ps.INK)
    spine.set_linewidth(1.6)

note = {"fontsize": 10, "color": ps.INK}
ax.text(x0 + 3.2 * sigma, y0, f"$(x_0,\\,y_0)=({x0:g},\\,{y0:g})$\n$\\sigma={sigma:g},\\ \\ A_0={a0:g}$", va="center",
        linespacing=1.6, **note)
ax.set_title(r"$\psi=\zeta=0,\ \ \partial c/\partial n=0\ \ \mathrm{on}\ \partial\Omega$", fontsize=10, pad=6)
ax.text(0.03, 0.96, "(a)", transform=ax.transAxes, va="top", **note)
ax.set_xlim(0, lx)
ax.set_ylim(0, ly)
ax.set_aspect("equal")
ax.set_xticks([0, lx / 2, lx])
ax.set_yticks([0, ly / 2, ly])
ax.xaxis.set_major_formatter("{x:g}")
ax.yaxis.set_major_formatter("{x:g}")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")

yy = np.linspace(0, ly, 200)
side.axvline(0, color=ps.GRID, linewidth=0.8)
side.plot(np.sin(2 * np.pi * yy / ly), yy, color=ps.BLUE, linewidth=1.4)
side.set_xlim(-1.25, 1.25)
side.set_xticks([-1, 0, 1])
side.set_xlabel("$F/F_0$")
side.set_title(r"$F=F_0\sin(2\pi y/L_y)$", fontsize=10, pad=6)
side.text(0.08, 0.96, "(b)", transform=side.transAxes, va="top", **note)
side.set_box_aspect(3.0)
ps.save(fig, "domain")
