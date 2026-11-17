#!/usr/bin/env python3
"""Report figure: the flow at release and the pigment at three times after it, from a snapshot folder.

Example:
    python3 scripts/plot_flow_evolution.py outputs/flow256 --steps 0 16000 64000
"""
import argparse
from pathlib import Path

import numpy as np

import fieldviz
import plot_style as ps

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("folder", help="snapshot folder (meta.json, index.csv, field_ and psi_ files)")
parser.add_argument("--steps", type=int, nargs=3, default=None, help="three snapshot steps (default: first, 1/4, last)")
parser.add_argument("--gamma", type=float, default=0.5, help="power of the colour scale (default 0.5, square root)")
parser.add_argument("--name", default="flow_evolution", help="file name in results/figures")
args = parser.parse_args()

folder = Path(args.folder)
meta = fieldviz.load_metadata(folder / "meta.json")
index = {row["step"]: row for row in fieldviz.read_snapshot_index(folder)}
steps = args.steps or [min(index), sorted(index)[len(index) // 4], max(index)]
shots = [index[s] for s in steps]
t0 = shots[0]["time"]
lx, ly = meta["lx"], meta["ly"]

ps.apply()
ps.plt.rcParams.update({"axes.grid": False, "axes.spines.top": True, "axes.spines.right": True, "mathtext.fontset": "cm"})
fig, axes = ps.plt.subplots(1, 4, figsize=(11.5, 3.2), layout="constrained", sharey=True)


def streamlines(ax, psi, alpha=1.0, width=0.6, arrows=False):
    """Thin grey contours of psi; with arrows, a few arrowheads show the direction of the flow."""
    up = psi[::-1]
    x = np.linspace(0, lx, up.shape[1])
    y = np.linspace(0, ly, up.shape[0])
    span = np.max(np.abs(up))
    levels = np.linspace(-span, span, 14)[1:-1]
    lines = ax.contour(x, y, up, levels=levels, colors=ps.MUTED, linewidths=width, linestyles="solid", alpha=alpha)
    if not arrows:
        return
    for level, segments in zip(lines.levels, lines.allsegs):
        if abs(abs(level) - 0.54 * span) > 0.08 * span:  # one contour in the middle of each gyre
            continue
        for path in segments:
            area = 0.5 * np.sum(path[:-1, 0] * path[1:, 1] - path[1:, 0] * path[:-1, 1])
            if (area > 0) != (level > 0):  # psi > 0 turns counter-clockwise
                path = path[::-1]
            for k in np.linspace(0, len(path) - 2, 4, endpoint=False).astype(int) + len(path) // 8:
                ax.annotate("", path[k + 1], path[k], arrowprops={"arrowstyle": "-|>", "color": ps.MUTED,
                                                                  "mutation_scale": 9, "linewidth": 0})


psi0 = fieldviz.load_field(shots[0]["psi"])
streamlines(axes[0], psi0, arrows=True)
axes[0].add_patch(ps.plt.Circle((meta["blob_x"], meta["blob_y"]), 2 * meta["blob_width"], facecolor=ps.BLUE,
                                alpha=0.25, edgecolor="none"))
axes[0].set_title("(a) Flow at release", fontsize=11)

peak = float(fieldviz.load_field(shots[0]["pigment"]).max())
norm, cmap = fieldviz.make_norm(0.0, peak, args.gamma)
letters = "bcd"
for ax, shot, letter in zip(axes[1:], shots, letters):
    image = fieldviz.draw_field(ax, fieldviz.load_field(shot["pigment"]), meta, norm, cmap)
    streamlines(ax, fieldviz.load_field(shot["psi"]), alpha=0.3, width=0.5)
    ax.set_title(f"({letter}) $\\tau = {shot['time'] - t0:.2f}$".replace("0.00", "0"), fontsize=11)

for ax in axes:
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    ax.set_aspect("equal")
    ax.set_xticks([0, lx / 2, lx])
    ax.set_yticks([0, ly / 2, ly])
    ax.xaxis.set_major_formatter("{x:g}")
    ax.yaxis.set_major_formatter("{x:g}")
    ax.set_xlabel("$x$")
axes[0].set_ylabel("$y$")
bar = fig.colorbar(image, ax=axes[1:], shrink=0.78, pad=0.02, ticks=[0, 0.05, 0.2, 0.4, 0.6, 0.8, 1.0], format="%g")
bar.set_label("pigment concentration $c$")
bar.outline.set_edgecolor(ps.MUTED)
ps.save(fig, args.name)
for shot in shots:
    c = fieldviz.load_field(shot["pigment"])
    print(f"step {shot['step']}: tau = {shot['time'] - t0:.3f}, max c = {c.max():.3f}, mass = {shot['mass']:.6f}")
