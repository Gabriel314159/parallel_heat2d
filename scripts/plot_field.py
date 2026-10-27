#!/usr/bin/env python3
"""Plot one scalar field written by the solver, optionally with the flow on top.

FIELD is a final field CSV (with its .json summary next to it) or one
snapshot .npy file from a snapshot folder (meta.json next to it).

Examples:
    ./build/heat_serial --output outputs/final.csv
    python3 scripts/plot_field.py outputs/final.csv --velocity streamlines

    python3 scripts/plot_field.py outputs/snapshots/field_005000.npy -o frame.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import fieldviz  # noqa: E402


def snapshot_time(field_path: Path, meta: dict) -> tuple[float, int | None, float | None]:
    """Time, step and mass of a snapshot file (from index.csv), or the final time of a CSV run."""
    index = field_path.parent / "index.csv"
    if field_path.suffix == ".npy" and index.is_file():
        for row in fieldviz.read_snapshot_index(field_path.parent):
            if row["file"].name == field_path.name:
                return row["time"], row["step"], row["mass"]
    return meta.get("t_end", 0.0), meta.get("steps"), meta.get("mass_final")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("field", help="field CSV or snapshot .npy")
    parser.add_argument("-o", "--output", help="image file (default: FIELD with .png)")
    parser.add_argument("--velocity", choices=["none", "streamlines", "arrows"], default="none",
                        help="draw the double-gyre flow at the field's time (default none)")
    parser.add_argument("--vmin", type=float, help="lower end of the colour scale (default: data minimum or 0)")
    parser.add_argument("--vmax", type=float, help="upper end of the colour scale (default: data maximum)")
    parser.add_argument("--gamma", type=float, default=1.0,
                        help="power colour scale for non-negative fields; <1 brightens faint pigment (default 1)")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    field_path = Path(args.field)
    field = fieldviz.load_field(field_path)
    meta = fieldviz.load_metadata(field_path)
    t, step, mass = snapshot_time(field_path, meta)

    vmin = args.vmin if args.vmin is not None else min(0.0, float(field.min()))
    vmax = args.vmax if args.vmax is not None else float(field.max())
    norm, cmap = fieldviz.make_norm(vmin, vmax, args.gamma)

    aspect = meta["ly"] / meta["lx"]
    fig, ax = plt.subplots(figsize=(6.4, 0.6 + 5.2 * aspect), facecolor=fieldviz.SURFACE, layout="constrained")
    image = fieldviz.draw_field(ax, field, meta, norm, cmap)
    fieldviz.draw_velocity(ax, meta, t, args.velocity)
    ax.set_xlim(0, meta["lx"])
    ax.set_ylim(0, meta["ly"])
    ax.set_aspect("equal")
    fieldviz.style_axes(ax)
    label = "concentration c" + (f"  (power scale, gamma = {args.gamma:g})" if args.gamma != 1.0 else "")
    fieldviz.style_colorbar(fig.colorbar(image, ax=ax, shrink=0.9), label)
    ax.set_title(fieldviz.title_text(t, step, mass), color=fieldviz.INK, fontsize=11, loc="left")
    fig.suptitle(fieldviz.subtitle_text(meta), color=fieldviz.INK_SECONDARY, fontsize=8.5, x=0.02, ha="left")

    output = Path(args.output) if args.output else field_path.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi, facecolor=fieldviz.SURFACE)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
