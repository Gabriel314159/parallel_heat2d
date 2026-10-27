#!/usr/bin/env python3
"""Run explicit experiment configurations; retain raw repetitions and provenance."""
import argparse
import copy
import csv
import datetime as dt
import itertools
import os
import sys
from common import (METRIC_COLUMNS, add_runtime_options, environment_record, load_experiment,
                    merge, numerical_parameters, resolve, run_configuration, save_json)


def configurations(experiment, base):
    kind = experiment['kind']
    if kind == 'correctness':
        raise ValueError('Use check_correctness.py for the correctness manifest')
    result = []
    if experiment.get('include_serial_reference'):
        result.append(('serial_reference', merge(base, {'execution': {'backend': 'serial', 'mpi_ranks': 1, 'omp_threads': 1}})))
    if kind == 'strategies':
        combinations = itertools.product(experiment['resources'], experiment['decompositions'], experiment['exchanges'])
    elif kind in ('strong_scaling', 'weak_scaling'):
        combinations = [(r, experiment['decomposition'], experiment['exchange']) for r in experiment['resources']]
    else:
        raise ValueError('Unknown experiment kind: ' + kind)
    for resource, decomposition, exchange in combinations:
        cfg = merge(base, {'execution': {**resource, 'backend': 'parallel', 'decomposition': decomposition, 'exchange': exchange}})
        if kind == 'strategies' and resource['mpi_ranks']*resource['omp_threads'] != experiment['total_cores']:
            raise ValueError('A strategy resource layout does not match total_cores')
        if kind == 'weak_scaling':
            py, px = resource['process_grid']
            p = cfg['problem']
            p['nx'], p['ny'] = experiment['local_nx']*px, experiment['local_ny']*py
            p['lx'] = (p['nx']+1)*experiment['spacing']['dx']
            p['ly'] = (p['ny']+1)*experiment['spacing']['dy']
        numerical_parameters(cfg['problem'])
        e = cfg['execution']
        result.append((f'p{e["mpi_ranks"]}_t{e["omp_threads"]}_{decomposition}_{exchange}', cfg))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='experiments/smoke_scaling.json')
    parser.add_argument('--repetitions', type=int)
    parser.add_argument('--warmup-runs', type=int)
    parser.add_argument('--allocated-cores', type=int, help='Explicit total allocated CPUs, including multi-node allocations')
    parser.add_argument('--dry-run', action='store_true', help='Write the resolved plan without running solvers')
    parser.add_argument('--keep-going', action='store_true')
    add_runtime_options(parser)
    args = parser.parse_args()
    path, experiment, base = load_experiment(args.config)
    runs = configurations(experiment, base)
    repetitions = args.repetitions if args.repetitions is not None else experiment['repetitions']
    warmups = args.warmup_runs if args.warmup_runs is not None else experiment['warmup_runs']
    if repetitions < 1 or warmups < 0:
        raise ValueError('repetitions must be positive; warmup-runs must be nonnegative')
    session = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = resolve(base['output']['directory'])/session
    metadata = environment_record(args.build_dir)
    metadata.update(config_file=str(path), session_id=session, launcher=args.launcher, mpi_args=args.mpi_arg,
                    repetitions=repetitions, warmup_runs=warmups, resolved_runs=[{'name': n, **c} for n, c in runs],
                    mpi_environment={k: os.environ[k] for k in ('MPIR_CVAR_CH4_NETMOD', 'MPIR_CVAR_CH4_CMA_ENABLE', 'FI_PROVIDER', 'UCX_TLS') if k in os.environ})
    save_json(out/'session.json', metadata)
    if args.dry_run:
        print(f'{len(runs)} configurations resolved. No solvers run. Plan: {out / "session.json"}')
        return
    available = args.allocated_cores if args.allocated_cores is not None else (len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1)
    if available < 1:
        raise ValueError('allocated-cores must be positive')
    for _, cfg in runs:
        e = cfg['execution']
        if e['mpi_ranks']*e['omp_threads'] > available:
            raise ValueError(f'Requested {e["mpi_ranks"]*e["omp_threads"]} cores but only {available} are available/declared')
    csv_path = resolve(experiment['metrics_csv'])
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    existed = csv_path.exists() and csv_path.stat().st_size > 0
    if existed:
        with csv_path.open(newline='') as f:
            if next(csv.reader(f)) != METRIC_COLUMNS:
                raise ValueError('Existing metrics CSV has a different schema; choose a new result path')
    failures = 0
    with csv_path.open('a', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=METRIC_COLUMNS, extrasaction='ignore')
        if not existed:
            writer.writeheader()
        for name, cfg in runs:
            for index in range(warmups+repetitions):
                phase = 'warmup' if index < warmups else 'measured'
                repetition = index if phase == 'warmup' else index-warmups
                row = {'experiment': path.stem, 'experiment_kind': experiment['kind'], 'session_id': session,
                       'run_id': name, 'phase': phase, 'repetition': repetition, 'status': 'ok',
                       'git_commit': metadata['git_commit'], 'timestamp_utc': dt.datetime.now(dt.timezone.utc).isoformat()}
                e, p = cfg['execution'], cfg['problem']
                row.update({k: p[k] for k in ('nx', 'ny', 'lx', 'ly', 'steps', 'alpha', 'case')})
                row.update({k: e[k] for k in ('backend', 'mpi_ranks', 'omp_threads', 'decomposition', 'exchange')})
                row['total_cores'] = e['mpi_ranks']*e['omp_threads']
                prefix = out/f'{name}_{phase}_{repetition}'
                field = prefix.with_suffix('.csv') if cfg['output'].get('write_field') else None
                try:
                    result = run_configuration(cfg, args, field, prefix)
                    if result['build_type'] not in ('Release', 'RelWithDebInfo'):
                        raise RuntimeError('Use a Release or RelWithDebInfo build for benchmarks')
                    row.update(result)
                    print(f'{name}: {phase} {repetition}, {result["elapsed_s"]:.6g} s')
                except Exception as error:
                    row.update(status='error', error=str(error))
                    failures += 1
                    writer.writerow(row); f.flush()
                    if not args.keep_going:
                        raise
                    continue
                writer.writerow(row); f.flush()
    print('Raw measurements:', csv_path)
    print('Resolved configuration and logs:', out)
    if failures:
        raise RuntimeError(f'{failures} solver invocations failed; failure records were retained')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('BENCHMARK FAILED:', error, file=sys.stderr)
        sys.exit(1)
