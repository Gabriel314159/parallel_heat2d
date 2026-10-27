"""Shared helpers for plot_field.py and animate_field.py (needs numpy and matplotlib).

Reads the solver's output formats (src/common/io.hpp):
  - a final field CSV plus its .json summary next to it,
  - a snapshot folder: meta.json, index.csv and field_NNNNNN.npy files.

Colour: a single-hue sequential scale (near-white = no pigment, dark blue =
concentrated) for non-negative fields, and a blue/red diverging scale with a
neutral grey at zero for signed fields such as cosine_decay.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, PowerNorm, TwoSlopeNorm

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
AXIS = "#c3c2b7"

SEQUENTIAL = LinearSegmentedColormap.from_list(
    "pigment", [SURFACE, "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
DIVERGING = LinearSegmentedColormap.from_list(
    "signed", ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f1a3a2", "#e34948", "#9c2b2b"])
STREAMLINE_COLOR = INK_SECONDARY


def load_field(path: str | Path) -> np.ndarray:
    path = Path(path)
    if path.suffix == ".npy":
        return np.load(path)
    return np.loadtxt(path, delimiter=",", ndmin=2)


def load_metadata(field_path: str | Path) -> dict:
    """The run summary next to a field CSV, or meta.json of a snapshot folder."""
    field_path = Path(field_path)
    for candidate in (field_path.with_suffix(".json"), field_path.parent / "meta.json"):
        if candidate.is_file():
            return json.loads(candidate.read_text())
    raise FileNotFoundError(f"no metadata for {field_path}: expected {field_path.with_suffix('.json').name} "
                            "or meta.json next to it")


def read_snapshot_index(folder: str | Path) -> list[dict]:
    with open(Path(folder) / "index.csv", newline="") as handle:
        return [{"step": int(row["step"]), "time": float(row["time"]), "mass": float(row["mass"]),
                 "file": Path(folder) / row["file"]} for row in csv.DictReader(handle)]


def make_norm(vmin: float, vmax: float, gamma: float = 1.0):
    """Diverging around 0 for signed data, otherwise sequential (optionally a power scale)."""
    span = max(abs(vmin), abs(vmax))
    if vmin < -1e-12 * span and vmax > 1e-12 * span:
        return TwoSlopeNorm(vcenter=0.0, vmin=-span, vmax=span), DIVERGING
    if vmax <= vmin:
        vmax = vmin + 1.0
    if gamma != 1.0:
        return PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax), SEQUENTIAL
    return Normalize(vmin=vmin, vmax=vmax), SEQUENTIAL


def velocity_on_grid(meta: dict, t: float, points: int = 40):
    """Analytic double-gyre velocity on an evenly spaced grid (y increasing), for overlays.

    psi = A sin(pi f(X, t)) sin(pi Y),  X = 2x/lx, Y = y/ly,  f = X + eps sin(wt)(X^2 - 2X)
    u = dpsi/dy,  v = -dpsi/dx   (docs/design.md, "Model").
    """
    lx, ly = meta["lx"], meta["ly"]
    amplitude, eps = meta["flow_amplitude"], meta["flow_epsilon"]
    omega = 2 * np.pi / meta["flow_period"]
    nx = max(8, int(round(points * lx / max(lx, ly))))
    ny = max(8, int(round(points * ly / max(lx, ly))))
    x = np.linspace(0, lx, nx)
    y = np.linspace(0, ly, ny)
    xx, yy = np.meshgrid(x, y)
    X = 2 * xx / lx
    Y = yy / ly
    s = eps * np.sin(omega * t)
    f = X + s * (X * X - 2 * X)
    df_dX = 1 + s * (2 * X - 2)
    u = amplitude * np.sin(np.pi * f) * (np.pi / ly) * np.cos(np.pi * Y)
    v = -amplitude * np.pi * np.cos(np.pi * f) * df_dX * (2 / lx) * np.sin(np.pi * Y)
    return x, y, u, v


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for spine in ax.spines.values():
        spine.set_color(AXIS)
        spine.set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=9)
    ax.set_xlabel("x", color=INK_SECONDARY)
    ax.set_ylabel("y", color=INK_SECONDARY)


def style_colorbar(colorbar, label: str) -> None:
    colorbar.outline.set_edgecolor(AXIS)
    colorbar.outline.set_linewidth(0.8)
    colorbar.ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=9)
    colorbar.set_label(label, color=INK_SECONDARY)


def draw_field(ax, field: np.ndarray, meta: dict, norm, cmap):
    """The field as cells on its physical extent; row 0 of the array is the top (y = ly)."""
    return ax.imshow(field, extent=(0, meta["lx"], 0, meta["ly"]), origin="upper", cmap=cmap, norm=norm,
                     interpolation="nearest")


def draw_velocity(ax, meta: dict, t: float, style: str):
    """Streamlines or arrows of the flow at time t.

    Returns every artist the overlay added, so an animation can remove them
    before drawing the next frame (streamplot adds its arrows as separate patches).
    """
    if style == "none" or meta.get("flow_amplitude", 0) == 0:
        return []
    before = set(ax.get_children())
    x, y, u, v = velocity_on_grid(meta, t, points=40 if style == "streamlines" else 20)
    if style == "streamlines":
        ax.streamplot(x, y, u, v, color=STREAMLINE_COLOR, linewidth=0.5, density=0.7, arrowsize=0.6,
                      broken_streamlines=False)
    elif style == "arrows":
        speed = float(np.max(np.hypot(u, v))) or 1.0
        ax.quiver(x, y, u, v, color=STREAMLINE_COLOR, width=0.0025, scale_units="xy", angles="xy",
                  scale=speed / (0.6 * (x[1] - x[0])))
    else:
        raise ValueError(f"unknown velocity style: {style}")
    added = [artist for artist in ax.get_children() if artist not in before]
    for artist in added:
        artist.set_alpha(0.3)
    return added


def title_text(t: float, step: int | None = None, mass: float | None = None) -> str:
    parts = [f"t = {t:.3f}"]
    if step is not None:
        parts.append(f"step {step}")
    if mass is not None:
        parts.append(f"total mass {mass:.6g}")
    return "   ".join(parts)


def subtitle_text(meta: dict) -> str:
    """Two short lines: the run, then the physics parameters."""
    grid = f"{meta['nx']} x {meta['ny']} cells"
    coarsen = meta.get("snapshot_coarsen", 1)
    if coarsen > 1:
        grid += f", shown as {coarsen} x {coarsen} block averages"
    return (f"{meta['case']}, {grid}\n"
            f"kappa = {meta['kappa']:g}, double gyre A = {meta['flow_amplitude']:g}, "
            f"eps = {meta['flow_epsilon']:g}, period = {meta['flow_period']:g}")
