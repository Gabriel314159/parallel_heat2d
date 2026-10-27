"""Shared configuration and subprocess helpers. All solver computation is C++."""
from __future__ import annotations
import copy
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
METRIC_COLUMNS = [
    'experiment', 'experiment_kind', 'session_id', 'run_id', 'repetition', 'phase', 'status',
    'backend', 'nx', 'ny', 'lx', 'ly', 'steps', 'alpha', 'dt', 'case',
    'mpi_ranks', 'omp_threads', 'total_cores', 'nodes', 'decomposition', 'exchange',
    'py', 'px', 'elapsed_s', 'hostname', 'compiler', 'mpi_library', 'build_type',
    'git_commit', 'timestamp_utc', 'error'
]


def resolve(path):
    p = Path(path).expanduser()
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save_json(path, value):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def merge(base, override):
    result = copy.deepcopy(base)
    for key, value in override.items():
        result[key] = merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else copy.deepcopy(value)
    return result


def load_experiment(path):
    p = resolve(path)
    experiment = read_json(p)
    if experiment.get('schema_version') != 1:
        raise ValueError('Unsupported experiment schema_version')
    defaults = read_json((p.parent / experiment['defaults']).resolve())
    base = merge(defaults, {k: experiment[k] for k in ('problem', 'execution', 'output') if k in experiment})
    return p, experiment, base


def numerical_parameters(problem):
    for key in ('nx', 'ny'):
        if type(problem[key]) is not int or problem[key] < 1:
            raise ValueError(f'{key} must be a positive integer')
    if type(problem['steps']) is not int or problem['steps'] < 0:
        raise ValueError('steps must be a nonnegative integer')
    for key in ('lx', 'ly', 'alpha'):
        if not math.isfinite(problem[key]) or problem[key] <= 0:
            raise ValueError(f'{key} must be finite and positive')
    if not 0 < problem['cfl'] <= 1:
        raise ValueError('cfl must be in (0,1]')
    dx, dy = problem['lx'] / (problem['nx'] + 1), problem['ly'] / (problem['ny'] + 1)
    step = problem.get('dt')
    if step is None:
        step = problem['cfl'] / (2 * problem['alpha'] * (1 / dx**2 + 1 / dy**2))
    if not math.isfinite(step) or step <= 0 or problem['alpha'] * step * (1 / dx**2 + 1 / dy**2) > .5 + 1e-14:
        raise ValueError('Unstable or invalid time step')
    return dx, dy, step


def problem_args(problem):
    numerical_parameters(problem)
    names = {'nx': '--nx', 'ny': '--ny', 'steps': '--steps', 'lx': '--lx', 'ly': '--ly',
             'alpha': '--alpha', 'cfl': '--cfl', 'dt': '--dt', 'case': '--case',
             'initial_temperature': '--initial-temperature', 'amplitude': '--amplitude'}
    args = []
    for key, option in names.items():
        if problem.get(key) is not None:
            args += [option, str(problem[key])]
    for key, value in problem.get('boundaries', {}).items():
        if key not in ('top', 'bottom', 'left', 'right'):
            raise ValueError('Unknown physical boundary: ' + key)
        args += ['--' + key, str(value)]
    return args


def add_runtime_options(parser):
    parser.add_argument('--build-dir', default='build')
    parser.add_argument('--launcher', default='mpiexec', help='MPI launcher executable')
    parser.add_argument('--mpi-arg', action='append', default=[], help='Repeat as --mpi-arg=VALUE for launcher options')
    parser.add_argument('--timeout', type=float, default=180, help='Timeout in seconds per solver invocation')


def run_configuration(config, options, field=None, log_prefix=None):
    problem, execution = config['problem'], config['execution']
    backend = execution['backend']
    ranks, threads = execution['mpi_ranks'], execution['omp_threads']
    if any(type(x) is not int or x < 1 for x in (ranks, threads)):
        raise ValueError('ranks and threads must be positive integers')
    if backend not in ('serial', 'parallel'):
        raise ValueError('backend must be serial or parallel')
    if backend == 'serial' and (ranks != 1 or threads != 1):
        raise ValueError('Serial runs require one rank and one thread')
    if config.get('output', {}).get('snapshot_interval', 0) != 0:
        raise ValueError('Intermediate snapshots are not implemented; set snapshot_interval to 0')
    binary = resolve(options.build_dir) / ('heat_serial' if backend == 'serial' else 'heat_parallel')
    if not binary.is_file():
        raise FileNotFoundError(f'Build the solver first: {binary}')
    args = [str(binary)] + problem_args(problem)
    if field is not None:
        args += ['--output', str(Path(field).resolve())]
    if backend == 'parallel':
        launcher = shutil.which(options.launcher)
        if launcher is None:
            raise FileNotFoundError('MPI launcher not found: ' + options.launcher)
        args += ['--decomp', execution['decomposition'], '--exchange', execution['exchange']]
        if execution.get('process_grid') is not None:
            py, px = execution['process_grid']
            if py * px != ranks:
                raise ValueError('Process grid does not match ranks')
            args += ['--py', str(py), '--px', str(px)]
        launcher_args = [token.replace('{threads}', str(threads)).replace('{ranks}', str(ranks)) for token in options.mpi_arg]
        args = [launcher, *launcher_args, '-n', str(ranks), *args]
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(threads)
    env['OMP_DYNAMIC'] = 'FALSE'
    completed = subprocess.run(args, cwd=ROOT, env=env, capture_output=True, text=True, timeout=options.timeout)
    if log_prefix is not None:
        p = Path(log_prefix)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.with_suffix('.stdout').write_text(completed.stdout, encoding='utf-8')
        p.with_suffix('.stderr').write_text(completed.stderr, encoding='utf-8')
    if completed.returncode:
        raise RuntimeError(f'Solver exited {completed.returncode}: {shlex.join(args)}\n{completed.stderr[-4000:]}')
    records = []
    for line in completed.stdout.splitlines():
        try:
            value = json.loads(line)
            if isinstance(value, dict) and 'elapsed_s' in value:
                records.append(value)
        except json.JSONDecodeError:
            pass
    if len(records) != 1:
        raise RuntimeError('Expected exactly one solver JSON result')
    result = records[0]
    if result['mpi_ranks'] != ranks or result['omp_threads'] != threads:
        raise RuntimeError('Observed rank/thread counts differ from the requested configuration')
    if not math.isfinite(result['elapsed_s']) or result['elapsed_s'] < 0:
        raise RuntimeError('Invalid solver time')
    return result


def environment_record(build_dir):
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True) if shutil.which('git') else None
    result = {'platform': platform.platform(), 'python': platform.python_version(),
              'logical_cpus': os.cpu_count(), 'git_commit': commit.stdout.strip() if commit and commit.returncode == 0 else 'uncommitted',
              'omp_places': os.environ.get('OMP_PLACES', ''), 'omp_proc_bind': os.environ.get('OMP_PROC_BIND', ''),
              'timestamp_utc': dt.datetime.now(dt.timezone.utc).isoformat()}
    if hasattr(os, 'sched_getaffinity'):
        result['cpu_affinity'] = sorted(os.sched_getaffinity(0))
    result['binary_sha256'] = {}
    for name in ('heat_serial', 'heat_parallel'):
        p = resolve(build_dir) / name
        if p.is_file():
            result['binary_sha256'][name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return result
