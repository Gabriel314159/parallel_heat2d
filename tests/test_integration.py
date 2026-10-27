#!/usr/bin/env python3
"""Extended MPI integration checks: all strategies, thin blocks and large messages."""
import argparse
import itertools
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from common import add_runtime_options, merge, read_json, resolve, run_configuration, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_runtime_options(parser)
    args = parser.parse_args()
    base = read_json(resolve('config/defaults.json'))
    out = resolve('outputs/extended_validation')
    out.mkdir(parents=True, exist_ok=True)
    checked = []
    for nx, ny, ranks in [(1,7,3), (7,1,3), (3,3,4), (7,5,6), (31,23,4)]:
        serial = merge(base, {'problem': {'nx':nx,'ny':ny,'steps':9}})
        run_configuration(serial, args, out/'serial.csv')
        expected = np.loadtxt(out/'serial.csv', delimiter=',')
        for decomposition, exchange, threads in itertools.product(('rows','2d'),('blocking','overlap'),(1,2)):
            if decomposition=='rows' and ranks>ny:
                continue
            cfg = merge(serial, {'execution': {'backend':'parallel','mpi_ranks':ranks,'omp_threads':threads,
                                               'decomposition':decomposition,'exchange':exchange}})
            run_configuration(cfg, args, out/'parallel.csv')
            actual = np.loadtxt(out/'parallel.csv', delimiter=',')
            np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-10)
            checked.append({'nx':nx,'ny':ny,'ranks':ranks,'threads':threads,'decomposition':decomposition,'exchange':exchange,
                            'max_abs_error':float(np.max(np.abs(actual-expected)))})
    # Rows/columns exceed typical eager-message limits; exercise rendezvous and derived datatypes.
    for nx, ny, grid, decomposition in [(8192,7,None,'rows'), (7,8192,[1,3],'2d')]:
        serial = merge(base, {'problem': {'nx':nx,'ny':ny,'steps':4}})
        run_configuration(serial,args,out/'large_serial.csv')
        expected = np.loadtxt(out/'large_serial.csv',delimiter=',')
        for exchange in ('blocking','overlap'):
            cfg=merge(serial,{'execution':{'backend':'parallel','mpi_ranks':3,'omp_threads':2,
                       'decomposition':decomposition,'exchange':exchange,'process_grid':grid}})
            run_configuration(cfg,args,out/'large_parallel.csv')
            actual=np.loadtxt(out/'large_parallel.csv',delimiter=',')
            np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-10)
            checked.append({'nx':nx,'ny':ny,'ranks':3,'threads':2,'decomposition':decomposition,'exchange':exchange,
                            'max_abs_error':float(np.max(np.abs(actual-expected)))})
    # All ranks should fail promptly for an impossible partition, rather than hang.
    command=[args.launcher,*args.mpi_arg,'-n','2',str(resolve(args.build_dir)/'heat_parallel'),
             '--nx','1','--ny','1','--decomp','2d']
    result=subprocess.run(command,capture_output=True,text=True,timeout=args.timeout)
    assert result.returncode!=0, 'Impossible process grid was accepted'
    save_json(out/'summary.json',{'passed':True,'cases':checked,'invalid_decomposition_rejected':True})
    print('PASS',len(checked),'extended serial/MPI comparisons and invalid-decomposition rejection')


if __name__=='__main__':
    main()
