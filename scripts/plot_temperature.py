#!/usr/bin/env python3
"""Plot a numeric field CSV including physical boundaries."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('field')
    parser.add_argument('--output', default='outputs/temperature.png')
    parser.add_argument('--lx', type=float, default=1)
    parser.add_argument('--ly', type=float, default=1)
    parser.add_argument('--vmin', type=float)
    parser.add_argument('--vmax', type=float)
    args = parser.parse_args()
    data = np.loadtxt(args.field, delimiter=',')
    if data.ndim != 2 or not np.isfinite(data).all() or args.lx <= 0 or args.ly <= 0:
        raise ValueError('Expected a finite 2D field and positive geometry')
    x = np.linspace(0, args.lx, data.shape[1])
    y = np.linspace(args.ly, 0, data.shape[0])
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    image = ax.pcolormesh(x, y, data, shading='nearest', cmap='inferno', vmin=args.vmin, vmax=args.vmax)
    ax.set(xlabel='x (model units)', ylabel='y (model units)', title='Temperature field', xlim=(0, args.lx), ylim=(0, args.ly))
    ax.set_aspect('equal')
    fig.colorbar(image, ax=ax, label='Temperature (°C)')
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(out)


if __name__ == '__main__':
    main()
