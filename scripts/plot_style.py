"""Shared rules for the report's curve plots (see results/figures/README.md)."""
import csv
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "results/raw"
FIGURES = REPO / "results/figures"

# Categorical slots in fixed order (never cycled); the same entity keeps its colour in every figure.
BLUE, ORANGE, AQUA, YELLOW, PURPLE = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cd1"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e3e2dd"
ROLE = {  # entity -> (colour, marker, linestyle)
    "serial": (BLUE, "o", "-"),          # the measured curve (speedup against heat_serial)
    "serial_omp": (ORANGE, "s", "-"),    # baseline heat_serial_omp
    "halo": (AQUA, "^", "-"),
    "sum": (YELLOW, "D", "-"),
    "rest": (PURPLE, "v", "-"),          # t_it minus the global sums: local sweeps and halo
    "ideal": (MUTED, None, "--"),
    "ideal_k": (MUTED, None, ":"),      # second reference: only the CG iterations grow
}


def apply():
    plt.rcParams.update({
        "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 10, "axes.titlesize": 13, "axes.titleweight": "normal", "axes.titlepad": 4, "axes.labelsize": 10,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
        "lines.linewidth": 1.8, "lines.markersize": 6, "lines.markeredgecolor": "white",
        "lines.markeredgewidth": 0.8,
        "legend.frameon": False, "legend.fontsize": 9,
    })


def load(*files, phase="measured"):
    """Rows of the ok, measured runs by run_id: {run_id: [row, ...]}."""
    runs = defaultdict(list)
    for name in files:
        for row in csv.DictReader((RAW / name).open(newline="")):
            if row["status"] == "ok" and row["phase"] == phase:
                runs[row["run_id"]].append(row)
    return runs


def stats(rows, key="elapsed_s"):
    """(median, min, max) of a column over the repetitions."""
    values = [float(r[key]) for r in rows]
    return st.median(values), min(values), max(values)


def per_iteration(rows, key="poisson_s"):
    """Median of key divided by the CG iterations of the run (steps x mean iterations per step)."""
    steps = float(rows[0]["steps"])
    return st.median(float(r[key]) / (steps * float(r["poisson_iterations_mean"])) for r in rows)


def curve(ax, x, y, role, label=None, lo=None, hi=None):
    colour, marker, style = ROLE[role]
    ax.plot(x, y, color=colour, marker=marker, linestyle=style, label=label, zorder=3)
    if lo is not None:
        ax.fill_between(x, lo, hi, color=colour, alpha=0.2, linewidth=0, zorder=2)


def ranks_axis(ax, ticks, label="MPI ranks (cores)"):
    """Log2 x axis with the plain integers as ticks."""
    ax.set_xscale("log", base=2)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(label)


def mark_nodes(ax, cores_per_node=64, label=True):
    """Dotted line where runs begin to span several nodes."""
    ax.axvline(cores_per_node, color=MUTED, linestyle=":", linewidth=1.2, zorder=1)
    if label:
        ax.annotate("> 1 node", (cores_per_node, 1), xycoords=("data", "axes fraction"),
                    xytext=(4, -12), textcoords="offset points", color=MUTED, fontsize=8)


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{ext}")
    plt.close(fig)
    print(f"results/figures/{name}.png and .pdf")
