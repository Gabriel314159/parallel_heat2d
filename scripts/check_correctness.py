#!/usr/bin/env python3
"""Numerical invariants, serial/MPI equivalence and independent sine refinement."""
import argparse
import copy
import math
import subprocess
import sys
import numpy as np
from common import (ROOT, add_runtime_options, environment_record, load_experiment, merge,
                    numerical_parameters, resolve, run_configuration, save_json)


def check_field(field, problem):
    if field.shape != (problem['ny'] + 2, problem['nx'] + 2) or not np.isfinite(field).all():
        raise AssertionError('Invalid field dimensions or nonfinite temperatures')
    b = problem['boundaries']
    for observed, expected in [(field[:, 0], b['left']), (field[:, -1], b['right']),
                               (field[0, 1:-1], b['top']), (field[-1, 1:-1], b['bottom'])]:
        np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=0)
    if problem['case'] == 'constant':
        np.testing.assert_allclose(field, problem['initial_temperature'], atol=1e-12, rtol=0)
    if problem['case'] == 'heated_top':
        lower = min(problem['initial_temperature'], *b.values())
        upper = max(problem['initial_temperature'], *b.values())
        if field.min() < lower - 1e-10 or field.max() > upper + 1e-10:
            raise AssertionError('Discrete maximum principle violated')


def serial_execution():
    return {'backend': 'serial', 'mpi_ranks': 1, 'omp_threads': 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='experiments/correctness.json')
    parser.add_argument('--out-dir', default='outputs/validation')
    parser.add_argument('--quick', action='store_true', help='Cap MPI cases at two ranks and one thread per rank')
    parser.add_argument('--serial-only', action='store_true')
    parser.add_argument('--skip-refinement', action='store_true')
    add_runtime_options(parser)
    args = parser.parse_args()
    _, experiment, base = load_experiment(args.config)
    out = resolve(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    checks = []
    # A hand-computed first step independently checks the shared stencil and indexing.
    hand = merge(base, {'problem': {'nx': 3, 'ny': 3, 'steps': 1}, 'execution': serial_execution()})
    run_configuration(hand, args, out/'first_step.csv')
    field = np.loadtxt(out/'first_step.csv', delimiter=',')
    np.testing.assert_allclose(field[1, 1:-1], 36.0, atol=1e-12)
    np.testing.assert_allclose(field[2:-1, 1:-1], 20.0, atol=1e-12)
    checks.append({'name': 'hand_computed_first_step', 'passed': True})
    for case in experiment['runs']:
        config = merge(base, case)
        name = case['name']
        serial = merge(config, {'execution': serial_execution()})
        serial_path = out/(name+'_serial.csv')
        serial_result = run_configuration(serial, args, serial_path, out/(name+'_serial'))
        reference = np.loadtxt(serial_path, delimiter=',')
        check_field(reference, serial['problem'])
        record = {'name': name, 'passed': True, 'nx': config['problem']['nx'], 'ny': config['problem']['ny']}
        if not args.serial_only:
            if args.quick:
                config['execution']['mpi_ranks'] = min(2, config['execution']['mpi_ranks'])
                config['execution']['omp_threads'] = 1
                config['execution']['process_grid'] = None
            parallel_path = out/(name+'_parallel.csv')
            result = run_configuration(config, args, parallel_path, out/(name+'_parallel'))
            actual = np.loadtxt(parallel_path, delimiter=',')
            check_field(actual, config['problem'])
            np.testing.assert_allclose(actual, reference,
                atol=experiment['checks']['serial_parallel_atol'], rtol=experiment['checks']['serial_parallel_rtol'])
            record.update(max_abs_error=float(np.max(np.abs(actual-reference))),
                          mpi_ranks=result['mpi_ranks'], omp_threads=result['omp_threads'])
        if config['problem']['case'] == 'sine_decay':
            p = config['problem']
            x = np.linspace(0, p['lx'], p['nx']+2)
            y = np.linspace(p['ly'], 0, p['ny']+2)
            exact = p.get('amplitude', 1) * np.sin(np.pi*y[:, None]/p['ly']) * np.sin(np.pi*x[None, :]/p['lx'])
            exact *= math.exp(-p['alpha']*math.pi**2*(1/p['lx']**2+1/p['ly']**2)*serial_result['final_time'])
            record['analytical_max_abs_error'] = float(np.max(np.abs(reference-exact)))
            np.testing.assert_allclose(record['analytical_max_abs_error'], serial_result['analytical_max_abs_error'], atol=1e-13)
        checks.append(record)
        print('PASS', name, 'max difference:', record.get('max_abs_error', 'serial only'))
    # Reject invalid numerical inputs, rather than silently running unstable problems.
    binary = resolve(args.build_dir)/'heat_serial'
    for options in [['--nx', '0'], ['--alpha', '0'], ['--alpha', 'nan'], ['--dt', '1'], ['--cfl', '1.5']]:
        completed = subprocess.run([str(binary), *options], capture_output=True, text=True, timeout=args.timeout)
        if completed.returncode == 0:
            raise AssertionError('Invalid parameters were accepted: ' + str(options))
    checks.append({'name': 'invalid_parameter_rejection', 'passed': True})
    if not args.skip_refinement:
        errors = []
        for n in (16, 32, 64):
            cfg = merge(base, {'problem': {'nx': n, 'ny': n, 'case': 'sine_decay', 'amplitude': 1.,
                        'boundaries': {'top': 0., 'bottom': 0., 'left': 0., 'right': 0.}}, 'execution': serial_execution()})
            _, _, safe_dt = numerical_parameters(cfg['problem'])
            cfg['problem']['steps'] = math.ceil(.01/safe_dt)
            cfg['problem']['dt'] = .01/cfg['problem']['steps']
            result = run_configuration(cfg, args)
            errors.append(result['analytical_max_abs_error'])
        if not (errors[1] < .6*errors[0] and errors[2] < .6*errors[1]):
            raise AssertionError('Sine-field errors did not converge under refinement: ' + str(errors))
        checks.append({'name': 'sine_refinement_at_t_0.01', 'passed': True, 'sizes': [16, 32, 64], 'errors': errors})
        print('PASS sine refinement:', errors)
    summary = {'passed': True, 'checks': checks, 'environment': environment_record(args.build_dir),
               'scope': 'Local numerical validation; not an HPC performance result.'}
    save_json(out/'summary.json', summary)
    print(f'All {len(checks)} checks passed. Summary: {out / "summary.json"}')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print('VALIDATION FAILED:', error, file=sys.stderr)
        sys.exit(1)
