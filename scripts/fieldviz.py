"""Shared helpers for plot_field.py and animate_field.py (needs numpy and matplotlib).

Reads the solver's output formats (src/common/io.hpp):
  - a final pigment field (CSV or .npy) plus its .json summary next to it,
  - a snapshot folder: meta.json, index.csv and, per snapshot, the pigment
    (field_NNNNNN.npy, cell values) and the streamfunction and vorticity
    (psi_NNNNNN.npy, vort_NNNNNN.npy, on the cell corners, walls included).

Colour: a single-hue sequential scale (near-white = no pigment, dark blue =
concentrated) for non-negative fields, and a blue/red diverging scale with a
neutral grey at zero for signed fields such as the vorticity.
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
OVERLAY_COLOR = INK_SECONDARY

FIELD_PREFIX = {"pigment": "field_", "vorticity": "vort_", "psi": "psi_"}
FIELD_LABEL = {"pigment": "concentration c", "vorticity": "vorticity ζ", "psi": "streamfunction ψ"}


def load_field(path: str | Path) -> np.ndarray:
    path = Path(path)
    if path.suffix == ".npy":
        return np.load(path)
    return np.loadtxt(path, delimiter=",", ndmin=2)


def load_metadata(field_path: str | Path) -> dict:
    """The run summary next to a field file, or meta.json of a snapshot folder."""
    field_path = Path(field_path)
    for candidate in (field_path.with_suffix(".json"), field_path.parent / "meta.json"):
        if candidate.is_file():
            return json.loads(candidate.read_text())
    raise FileNotFoundError(f"no metadata for {field_path}: expected {field_path.with_suffix('.json').name} "
                            "or meta.json next to it")


def read_snapshot_index(folder: str | Path) -> list[dict]:
    """One dict per snapshot: step, time, mass, energy, enstrophy and the paths of its files."""
    folder = Path(folder)
    rows = []
    with open(folder / "index.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            entry = {"step": int(row["step"]), "time": float(row["time"]), "mass": float(row["mass"])}
            for key in ("energy", "enstrophy"):
                if key in row:
                    entry[key] = float(row[key])
            entry["pigment"] = folder / row["file"]
            if "psi_file" in row:
                entry["psi"] = folder / row["psi_file"]
                entry["vorticity"] = folder / row["vort_file"]
            rows.append(entry)
    return rows


def companion(path: Path, kind: str) -> Path:
    """The snapshot file of another kind with the same step, e.g. psi_000100.npy for field_000100.npy."""
    for prefix in FIELD_PREFIX.values():
        if path.name.startswith(prefix):
            return path.with_name(FIELD_PREFIX[kind] + path.name[len(prefix):])
    raise ValueError(f"{path.name} is not a snapshot file")


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
    """Cell values (pigment) or corner values (psi, vorticity) on the basin; row 0 is the top (y = ly)."""
    return ax.imshow(field, extent=(0, meta["lx"], 0, meta["ly"]), origin="upper", cmap=cmap, norm=norm,
                     interpolation="nearest" if field.shape[0] > 64 else "bilinear")


def draw_flow(ax, psi: np.ndarray | None, meta: dict, style: str):
    """Streamlines (contours of psi) or arrows of the flow, from a psi snapshot on the corners.

    Returns every artist added, so that an animation can remove them before the next frame.
    """
    if style == "none" or psi is None or not np.any(psi):
        return []
    before = set(ax.get_children())
    ny, nx = psi.shape[0] - 1, psi.shape[1] - 1
    x = np.linspace(0, meta["lx"], nx + 1)
    y = np.linspace(0, meta["ly"], ny + 1)
    up = psi[::-1]  # rows from y = 0 upwards
    if style == "streamlines":
        span = np.max(np.abs(up))
        levels = np.linspace(-span, span, 22)[1:-1]
        ax.contour(x, y, up, levels=levels, colors=OVERLAY_COLOR, linewidths=0.6, linestyles="solid")
    elif style == "arrows":
        u = np.gradient(up, y, axis=0)    # u = dpsi/dy
        v = -np.gradient(up, x, axis=1)   # v = -dpsi/dx
        step = max(1, nx // 20)
        speed = float(np.max(np.hypot(u, v))) or 1.0
        ax.quiver(x[::step], y[::step], u[::step, ::step], v[::step, ::step], color=OVERLAY_COLOR,
                  width=0.0025, scale_units="xy", angles="xy", scale=speed / (0.9 * step * (x[1] - x[0])))
    else:
        raise ValueError(f"unknown flow style: {style}")
    added = [artist for artist in ax.get_children() if artist not in before]
    for artist in added:
        artist.set_alpha(0.35)
    return added


def subtitle_text(meta: dict) -> str:
    """Two short lines: the run, then the flow parameters."""
    grid = f"{meta['nx']} x {meta['ny']} cells"
    coarsen = meta.get("snapshot_coarsen", 1)
    if coarsen > 1:
        grid += f", shown {coarsen}x coarser"
    wind = f"wind {meta['wind']} F0 = {meta['wind_amplitude']:g}"
    if meta.get("wind") == "periodic":
        wind += f", alpha = {meta['wind_alpha']:g}, period = {meta['wind_period']:g}"
    return (f"{meta['case']}, {grid}\n"
            f"viscosity = {meta['viscosity']:g}, {wind}, kappa = {meta['kappa']:g}")
