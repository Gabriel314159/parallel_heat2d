#!/usr/bin/env python3
"""Compare a serial run with a many-rank run through the job's launcher, field by field."""
import argparse
from pathlib import Path
import sys
from runner_common import run_case, read_field


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build-dir", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--launcher", required=True)
    p.add_argument("--ranks", type=int, default=256)
    p.add_argument("--n", type=int, default=48)
    args = p.parse_args()
    out = Path(args.output).resolve()
    base = dict(nx=args.n, ny=args.n, nz=args.n, lx=1, ly=1, lz=1, steps=5, spinup=0, tol=1e-10)
    serial, parallel = out / "serial.csv", out / "parallel.csv"
    run_case(args.build_dir, args.launcher, out, "serial", dict(base, output=serial))
    r = run_case(args.build_dir, args.launcher, out, "parallel", dict(base, output=parallel), args.ranks, 1)
    a, b = read_field(serial), read_field(parallel)
    if a.keys() != b.keys():
        raise RuntimeError("the two runs do not hold the same cells")
    for field in ("c", "phi", "u_e", "v_n", "w_t"):
        delta = max(abs(a[key][field] - b[key][field]) for key in a)
        limit = 1e-8 if field == "phi" else 2e-10
        print("{} max difference {:.3e} (limit {:.0e})".format(field, delta, limit), flush=True)
        if not delta <= limit:
            raise RuntimeError("serial and {}-rank runs differ in {}".format(args.ranks, field))
    print("PASS: serial equals {} ranks on {} nodes, process grid {} x {} x {}".format(
        r["ranks"], r["nodes"], r["px"], r["py"], r["pz"]), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        sys.exit(1)
