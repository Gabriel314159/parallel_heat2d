#!/usr/bin/env python3
"""Plot measured medians and scaling, using one experiment session at a time."""
import argparse
import collections
import csv
from pathlib import Path
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv')
    parser.add_argument('--out-dir', default='outputs/performance_plots')
    parser.add_argument('--session', default='latest')
    args = parser.parse_args()
    with open(args.csv, newline='', encoding='utf-8') as f:
        rows = [r for r in csv.DictReader(f) if r['status'] == 'ok' and r['phase'] == 'measured']
    if not rows:
        raise ValueError('No successful measured records; a CSV header is not experimental data')
    session = max(r['session_id'] for r in rows) if args.session == 'latest' else args.session
    rows = [r for r in rows if r['session_id'] == session]
    if not rows:
        raise ValueError('Requested session not found')
    kinds = {r['experiment_kind'] for r in rows}
    if len(kinds) != 1:
        raise ValueError('Plot one experiment kind per session')
    kind = kinds.pop()
    # Guard against combining unlike discrete problems into a strong-scaling baseline.
    if kind != 'weak_scaling':
        signatures = {tuple(r[k] for k in ('nx','ny','lx','ly','steps','alpha','dt','case')) for r in rows}
        if len(signatures) != 1:
            raise ValueError('This session contains different numerical problems')
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r['run_id']].append(r)
    points = []
    for name, samples in groups.items():
        values = [float(r['elapsed_s']) for r in samples]
        if min(values) <= 0 or not np.isfinite(values).all():
            raise ValueError('Scaling requires positive finite timings')
        points.append({**samples[0], 'median': statistics.median(values), 'minimum': min(values), 'maximum': max(values)})
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if kind == 'strategies':
        strategies = sorted({(p['decomposition'], p['exchange']) for p in points})
        layouts = sorted({(int(p['mpi_ranks']), int(p['omp_threads'])) for p in points}, reverse=True)
        fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
        x = np.arange(len(layouts))
        width = .8/len(strategies)
        for index, strategy in enumerate(strategies):
            lookup = {(int(p['mpi_ranks']), int(p['omp_threads'])): p['median'] for p in points if (p['decomposition'], p['exchange']) == strategy}
            ax.bar(x-.4+width*(index+.5), [lookup.get(layout, np.nan) for layout in layouts], width=width, label='/'.join(strategy))
        ax.set_xticks(x, [f'{r} × {t}' for r,t in layouts])
        ax.set(xlabel='MPI ranks × threads per rank', ylabel='Median solver time (s)', title='Fixed-resource strategy comparison')
        ax.legend()
        fig.savefig(out/'strategies.png', dpi=180); plt.close(fig)
    else:
        series = collections.defaultdict(list)
        serial = [p for p in points if p['backend'] == 'serial']
        for p in points:
            if p['backend'] == 'parallel':
                series[(p['decomposition'], p['exchange'], p['omp_threads'])].append(p)
        figures = ['runtime']
        if kind == 'strong_scaling' and serial:
            if len(serial) != 1:
                raise ValueError('Expected exactly one matching serial reference')
            figures += ['speedup', 'efficiency']
        if kind == 'weak_scaling':
            figures += ['weak_efficiency']
        for metric in figures:
            fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
            for label, samples in sorted(series.items()):
                samples.sort(key=lambda p:int(p['total_cores']))
                x = [int(p['total_cores']) for p in samples]
                times = np.array([p['median'] for p in samples])
                if metric == 'runtime':
                    y = times
                    error = np.array([[p['median']-p['minimum'] for p in samples], [p['maximum']-p['median'] for p in samples]])
                    ax.errorbar(x,y,yerr=error,marker='o',capsize=3,label=f'{label[0]}/{label[1]}, {label[2]} thread(s)/rank')
                else:
                    reference = serial[0]['median'] if metric in ('speedup', 'efficiency') else times[0]
                    y = reference/times
                    if metric == 'efficiency': y /= np.array(x)
                    ax.plot(x,y,'o-',label=f'{label[0]}/{label[1]}, {label[2]} thread(s)/rank')
            if metric == 'runtime' and serial:
                ax.axhline(serial[0]['median'],color='gray',linestyle='--',label='Pure serial reference')
            ylabel = {'runtime':'Solver time (s), median and min–max', 'speedup':'Speedup vs pure serial',
                      'efficiency':'Strong-scaling efficiency', 'weak_efficiency':'Efficiency vs smallest same-thread MPI layout'}[metric]
            ax.set(xlabel='Total CPU cores (ranks × threads)',ylabel=ylabel,title=metric.replace('_',' ').title())
            ax.grid(alpha=.2); ax.legend(fontsize=8)
            fig.savefig(out/(metric+'.png'),dpi=180); plt.close(fig)
    print('Session:',session)
    print('Plots:',out)


if __name__ == '__main__':
    main()
