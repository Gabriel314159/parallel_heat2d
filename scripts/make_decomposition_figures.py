#!/usr/bin/env python3
"""Draw the two domain decompositions of docs/design.md as small figures.

Writes docs/diagrams/decomposition_rows.png and decomposition_blocks.png.
Needs NumPy-free Matplotlib only. Run from the repository root:
    python3 scripts/make_decomposition_figures.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

OUT = Path(__file__).resolve().parent.parent / "docs" / "diagrams"
COLORS = ["#cfe2f3", "#fce5cd", "#d9ead3", "#ead1dc", "#fff2cc", "#d0e0e3"]
INK = "#333333"


def split(n, p):
    """Sizes of n items over p parts: the first n % p parts get one extra."""
    q, r = divmod(n, p)
    return [q + 1 if k < r else q for k in range(p)]


def draw_blocks(ax, ny, nx, py, px):
    """Cells of an ny x nx grid (row 0 on top), coloured per rank."""
    rows, cols = split(ny, py), split(nx, px)
    i0 = 0
    rank = 0
    for a, h in enumerate(rows):
        j0 = 0
        for b, w in enumerate(cols):
            ax.add_patch(Rectangle((j0, ny - i0 - h), w, h, facecolor=COLORS[rank % 6],
                                   edgecolor=INK, linewidth=2.2))
            for i in range(i0, i0 + h):
                for j in range(j0, j0 + w):
                    ax.add_patch(Rectangle((j, ny - i - 1), 1, 1, fill=False,
                                           edgecolor=INK, linewidth=0.3, alpha=0.5))
            ax.text(j0 + w / 2, ny - i0 - h / 2, f"rank {rank}\n{h} x {w}",
                    ha="center", va="center", fontsize=10, color=INK)
            j0 += w
            rank += 1
        i0 += h
    return rows, cols


def arrow(ax, p, q):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<->", mutation_scale=12,
                                 color="#b00020", linewidth=1.6, zorder=5))


def finish(ax, ny, nx, title):
    ax.set_xlim(-0.6, nx + 0.6)
    ax.set_ylim(-0.6, ny + 0.6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=11, color=INK)


def rows_figure():
    ny, nx, p = 10, 8, 4
    fig, ax = plt.subplots(figsize=(4.2, 4.4))
    rows, _ = draw_blocks(ax, ny, nx, p, 1)
    edge = 0
    for h in rows[:-1]:
        edge += h
        y = ny - edge
        for x in (nx * 0.3, nx * 0.7):
            arrow(ax, (x, y + 0.45), (x, y - 0.45))
    finish(ax, ny, nx, "Rows: 10 x 8 cells, 4 ranks (3, 3, 2, 2 rows)")
    ax.text(nx / 2, -0.4, "red arrows: halo rows exchanged with the neighbour above and below",
            ha="center", fontsize=7.5, color="#b00020")
    return fig


def blocks_figure():
    ny, nx, py, px = 10, 8, 2, 3
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    rows, cols = draw_blocks(ax, ny, nx, py, px)
    y = ny - rows[0]
    x = 0
    for w in cols:
        arrow(ax, (x + w / 2, y + 0.45), (x + w / 2, y - 0.45))
        x += w
    x = 0
    for w in cols[:-1]:
        x += w
        for yy in (ny - rows[0] / 2 + 1.2, ny - rows[0] - rows[1] / 2 + 1.2):
            arrow(ax, (x - 0.45, yy), (x + 0.45, yy))
    finish(ax, ny, nx, "2D blocks: 10 x 8 cells, 2 x 3 process grid")
    ax.text(nx / 2, -0.4, "red arrows: halo rows and columns exchanged with each neighbour",
            ha="center", fontsize=7.5, color="#b00020")
    return fig


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for name, make in (("decomposition_rows", rows_figure), ("decomposition_blocks", blocks_figure)):
        fig = make()
        fig.savefig(OUT / f"{name}.png", dpi=170, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print("wrote", OUT / f"{name}.png")
