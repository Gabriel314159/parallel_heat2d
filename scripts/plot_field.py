#!/usr/bin/env python3
"""Plot one field written by the solver, optionally with the streamlines of the flow.

FIELD is a final pigment field (.csv or .npy, with its .json summary next to
it) or one snapshot file (field_, psi_ or vort_NNNNNN.npy, with meta.json next
to it). For snapshots, --field picks what to show and --flow draws the flow of
the same moment from its psi file.

Examples:
    python3 scripts/plot_field.py outputs/preview/field_004000.npy --flow streamlines
    python3 scripts/plot_field.py outputs/preview/field_004000.npy --field vorticity
    python3 scripts/plot_field.py outputs/final.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import fieldviz  # noqa: E402


def snapshot_info(path: Path, meta: dict) -> tuple[float, int | None]:
    """Time and step of a snapshot file (from index.csv), or the final time of a run summary."""
    index = path.parent / "index.csv"
    if path.suffix == ".npy" and index.is_file():
        for row in fieldviz.read_snapshot_index(path.parent):
            if path in (row["pigment"], row.get("psi"), row.get("vorticity")):
                return row["time"], row["step"]
    return meta.get("t_end", 0.0), meta.get("steps")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path", metavar="FIELD", help="final field (.csv/.npy) or a snapshot file")
    parser.add_argument("-o", "--output", help="image file (default: FIELD with .png)")
    parser.add_argument("--field", choices=["pigment", "vorticity", "psi"],
                        help="for a snapshot: what to show (default: the kind of FIELD)")
    parser.add_argument("--flow", choices=["none", "streamlines", "arrows"], default="none",
                        help="for a snapshot: draw the flow from its psi file (default none)")
    parser.add_argument("--vmin", type=float, help="lower end of the colour scale")
    parser.add_argument("--vmax", type=float, help="upper end of the colour scale")
    parser.add_argument("--gamma", type=float, default=1.0,
                        help="power colour scale for the pigment; <1 brightens faint pigment (default 1)")
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    path = Path(args.path)
    meta = fieldviz.load_metadata(path)
    is_snapshot = (path.parent / "index.csv").is_file() and path.suffix == ".npy"
    kind = args.field or next((k for k, p in fieldviz.FIELD_PREFIX.items() if path.name.startswith(p)), "pigment")
    if is_snapshot and args.field:
        path = fieldviz.companion(path, kind)
    field = fieldviz.load_field(path)
    psi = None
    if is_snapshot and args.flow != "none":
        psi = fieldviz.load_field(fieldviz.companion(path, "psi"))
    elif args.flow != "none":
        print("note: --flow needs a snapshot (psi is not saved with final fields); drawing no flow", file=sys.stderr)
    t, step = snapshot_info(path, meta)

    vmin = args.vmin if args.vmin is not None else float(min(0.0, field.min()))
    vmax = args.vmax if args.vmax is not None else float(field.max())
    norm, cmap = fieldviz.make_norm(vmin, vmax, args.gamma if kind == "pigment" else 1.0)

    aspect = meta["ly"] / meta["lx"]
    fig, ax = plt.subplots(figsize=(4.2 + 2.0 / max(aspect, 0.5), 0.9 + 4.6 * min(aspect, 2.0) / 1.0 * 0.75),
                           facecolor=fieldviz.SURFACE, layout="constrained")
    image = fieldviz.draw_field(ax, field, meta, norm, cmap)
    fieldviz.draw_flow(ax, psi, meta, args.flow)
    ax.set_xlim(0, meta["lx"])
    ax.set_ylim(0, meta["ly"])
    ax.set_aspect("equal")
    fieldviz.style_axes(ax)
    label = fieldviz.FIELD_LABEL[kind]
    if kind == "pigment" and args.gamma != 1.0:
        label += f"  (power scale, gamma = {args.gamma:g})"
    fieldviz.style_colorbar(fig.colorbar(image, ax=ax, shrink=0.9), label)
    title = f"t = {t:.2f}" + (f"   step {step}" if step is not None else "")
    ax.set_title(title, color=fieldviz.INK, fontsize=11, loc="left")
    fig.suptitle(fieldviz.subtitle_text(meta), color=fieldviz.INK_SECONDARY, fontsize=8.5, x=0.02, ha="left")

    output = Path(args.output) if args.output else path.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi, facecolor=fieldviz.SURFACE)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    np.seterr(all="ignore")
    sys.exit(main())
