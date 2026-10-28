#!/usr/bin/env python3
"""Animate a snapshot folder written with --snapshot-interval.

--field picks what to show: the pigment (default), the vorticity or the
streamfunction. The colour scale is fixed over the whole animation. For the
pigment, no frame leaves the range of the first (discrete maximum principle),
and the default power scale (gamma 0.5) keeps faint, spread-out pigment
visible; use --gamma 1 for a linear scale. --flow draws the streamlines of
each frame from its psi file.

Examples:
    ./build/heat_serial --snapshot-interval 500 --snapshot-dir outputs/preview
    python3 scripts/animate_field.py outputs/preview --flow streamlines
    python3 scripts/animate_field.py outputs/preview --field vorticity -o outputs/vorticity.mp4

The output format follows the extension: .gif (Pillow) or .mp4 (needs ffmpeg).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter  # noqa: E402

import fieldviz  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folder", help="snapshot folder (contains meta.json, index.csv and the .npy files)")
    parser.add_argument("-o", "--output", help="animation file, .gif or .mp4 (default: FOLDER_FIELD.gif)")
    parser.add_argument("--field", choices=["pigment", "vorticity", "psi"], default="pigment")
    parser.add_argument("--flow", choices=["none", "streamlines", "arrows"], default="none",
                        help="draw the flow of each frame (default none)")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--every", type=int, default=1, help="use every k-th snapshot (default 1)")
    parser.add_argument("--gamma", type=float, default=0.5,
                        help="power colour scale for the pigment; 1 = linear (default 0.5)")
    parser.add_argument("--vmax", type=float, help="upper end of the colour scale (default: from the frames)")
    parser.add_argument("--dpi", type=int, default=110)
    parser.add_argument("--frames-dir", help="also save every frame as PNG in this folder")
    args = parser.parse_args()

    folder = Path(args.folder)
    meta = json.loads((folder / "meta.json").read_text())
    frames = fieldviz.read_snapshot_index(folder)[:: max(1, args.every)]
    if not frames:
        print(f"no snapshots listed in {folder / 'index.csv'}", file=sys.stderr)
        return 1
    if args.field != "pigment" and "psi" not in frames[0]:
        print("this snapshot folder has no psi/vorticity files", file=sys.stderr)
        return 1

    first = fieldviz.load_field(frames[0][args.field])
    if args.field == "pigment":
        vmin = min(0.0, float(first.min()))
        vmax = args.vmax if args.vmax is not None else float(first.max())
        norm, cmap = fieldviz.make_norm(vmin, vmax, args.gamma)
    else:
        # A symmetric scale large enough for every frame (they start at rest).
        span = args.vmax if args.vmax is not None else max(
            float(np.max(np.abs(np.load(f[args.field], mmap_mode="r")))) for f in frames)
        norm, cmap = fieldviz.make_norm(-span or -1.0, span or 1.0)

    aspect = meta["ly"] / meta["lx"]
    fig, ax = plt.subplots(figsize=(4.2 + 2.0 / max(aspect, 0.5), 0.9 + 4.6 * min(aspect, 2.0) * 0.75),
                           facecolor=fieldviz.SURFACE, layout="constrained")
    image = fieldviz.draw_field(ax, first, meta, norm, cmap)
    ax.set_xlim(0, meta["lx"])
    ax.set_ylim(0, meta["ly"])
    ax.set_aspect("equal")
    fieldviz.style_axes(ax)
    label = fieldviz.FIELD_LABEL[args.field]
    if args.field == "pigment" and args.gamma != 1.0:
        label += f"  (power scale, gamma = {args.gamma:g})"
    fieldviz.style_colorbar(fig.colorbar(image, ax=ax, shrink=0.9), label)
    fig.suptitle(fieldviz.subtitle_text(meta), color=fieldviz.INK_SECONDARY, fontsize=8.5, x=0.02, ha="left")
    title = ax.set_title("", color=fieldviz.INK, fontsize=11, loc="left")
    overlay: list = []
    mass0 = frames[0]["mass"]

    frames_dir = Path(args.frames_dir) if args.frames_dir else None
    if frames_dir:
        frames_dir.mkdir(parents=True, exist_ok=True)

    def draw(k: int):
        frame = frames[k]
        image.set_data(fieldviz.load_field(frame[args.field]))
        for artist in overlay:
            artist.remove()
        psi = fieldviz.load_field(frame["psi"]) if args.flow != "none" and "psi" in frame else None
        overlay[:] = fieldviz.draw_flow(ax, psi, meta, args.flow)
        text = f"t = {frame['time']:.1f}   step {frame['step']}"
        if args.field == "pigment":
            change = abs(frame["mass"] - mass0) / abs(mass0) if mass0 else abs(frame["mass"] - mass0)
            text += f"   mass change {change:.1e}"
        elif "energy" in frame:
            text += f"   energy {frame['energy']:.3g}"
        title.set_text(text)
        if frames_dir:
            fig.savefig(frames_dir / f"frame_{k:04d}.png", dpi=args.dpi, facecolor=fieldviz.SURFACE)
        return [image, title, *overlay]

    output = Path(args.output) if args.output else folder.with_name(f"{folder.name}_{args.field}.gif")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == ".mp4":
        writer = FFMpegWriter(fps=args.fps, bitrate=2400)
    elif output.suffix == ".gif":
        writer = PillowWriter(fps=args.fps)
    else:
        print("output must end in .gif or .mp4", file=sys.stderr)
        return 1

    animation = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    animation.save(output, writer=writer, dpi=args.dpi, savefig_kwargs={"facecolor": fieldviz.SURFACE})
    print(f"wrote {output} ({len(frames)} frames, t = {frames[0]['time']:.1f} .. {frames[-1]['time']:.1f})")
    return 0


if __name__ == "__main__":
    np.seterr(all="ignore")
    sys.exit(main())
