#!/usr/bin/env python3
"""Animate a snapshot folder written with --snapshot-interval.

The colour scale is fixed over the whole animation, so fading of the pigment
by diffusion is visible. By the discrete maximum principle no frame exceeds
the range of the first one. The default power scale (gamma 0.5) keeps faint,
spread-out pigment visible; use --gamma 1 for a linear scale.

Examples:
    ./build/heat_serial --steps 20400 --snapshot-interval 204 --snapshot-dir outputs/preview
    python3 scripts/animate_field.py outputs/preview --velocity streamlines
    python3 scripts/animate_field.py outputs/preview -o outputs/preview.mp4 --fps 25

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
    parser.add_argument("folder", help="snapshot folder (contains meta.json, index.csv, field_*.npy)")
    parser.add_argument("-o", "--output", help="animation file, .gif or .mp4 (default: FOLDER.gif)")
    parser.add_argument("--velocity", choices=["none", "streamlines", "arrows"], default="none",
                        help="draw the double-gyre flow of each frame (default none)")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--every", type=int, default=1, help="use every k-th snapshot (default 1)")
    parser.add_argument("--gamma", type=float, default=0.5,
                        help="power colour scale for non-negative fields; 1 = linear (default 0.5)")
    parser.add_argument("--vmax", type=float, help="upper end of the colour scale (default: first frame maximum)")
    parser.add_argument("--dpi", type=int, default=110)
    parser.add_argument("--frames-dir", help="also save every frame as PNG in this folder")
    args = parser.parse_args()

    folder = Path(args.folder)
    meta = json.loads((folder / "meta.json").read_text())
    frames = fieldviz.read_snapshot_index(folder)[:: max(1, args.every)]
    if not frames:
        print(f"no snapshots listed in {folder / 'index.csv'}", file=sys.stderr)
        return 1

    first = fieldviz.load_field(frames[0]["file"])
    vmin = min(0.0, float(first.min()))
    vmax = args.vmax if args.vmax is not None else float(first.max())
    norm, cmap = fieldviz.make_norm(vmin, vmax, args.gamma)
    mass0 = frames[0]["mass"]

    aspect = meta["ly"] / meta["lx"]
    fig, ax = plt.subplots(figsize=(6.4, 0.6 + 5.2 * aspect), facecolor=fieldviz.SURFACE, layout="constrained")
    image = fieldviz.draw_field(ax, first, meta, norm, cmap)
    ax.set_xlim(0, meta["lx"])
    ax.set_ylim(0, meta["ly"])
    ax.set_aspect("equal")
    fieldviz.style_axes(ax)
    label = "concentration c" + (f"  (power scale, gamma = {args.gamma:g})" if args.gamma != 1.0 else "")
    fieldviz.style_colorbar(fig.colorbar(image, ax=ax, shrink=0.9), label)
    fig.suptitle(fieldviz.subtitle_text(meta), color=fieldviz.INK_SECONDARY, fontsize=8.5, x=0.02, ha="left")
    title = ax.set_title("", color=fieldviz.INK, fontsize=11, loc="left")
    overlay: list = []

    frames_dir = Path(args.frames_dir) if args.frames_dir else None
    if frames_dir:
        frames_dir.mkdir(parents=True, exist_ok=True)

    def draw(k: int):
        frame = frames[k]
        image.set_data(fieldviz.load_field(frame["file"]))
        for artist in overlay:
            artist.remove()
        overlay[:] = fieldviz.draw_velocity(ax, meta, frame["time"], args.velocity)
        change = abs(frame["mass"] - mass0) / abs(mass0) if mass0 else abs(frame["mass"] - mass0)
        title.set_text(f"t = {frame['time']:.2f}   step {frame['step']}   mass change {change:.1e}")
        if frames_dir:
            fig.savefig(frames_dir / f"frame_{k:04d}.png", dpi=args.dpi, facecolor=fieldviz.SURFACE)
        return [image, title, *overlay]

    output = Path(args.output) if args.output else folder.with_suffix(".gif")
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
    print(f"wrote {output} ({len(frames)} frames, t = {frames[0]['time']:.2f} .. {frames[-1]['time']:.2f})")
    return 0


if __name__ == "__main__":
    np.seterr(all="ignore")
    sys.exit(main())
