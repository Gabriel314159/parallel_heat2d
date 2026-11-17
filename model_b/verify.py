#!/usr/bin/env python3
"""Verify Neumann pressure, projection, pigment, and distributed full-field agreement."""
import argparse
import json
import math
from pathlib import Path
import sys
from runner_common import run_case, read_field, write_csv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build-dir", default="build")
    p.add_argument("--output", default="validation")
    p.add_argument("--launcher", default="mpirun --bind-to none")
    p.add_argument("--serial-only", action="store_true")
    p.add_argument("--timeout", type=int, default=600)
    args = p.parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out/"status.json").write_text(json.dumps(dict(status="RUNNING"))+"\n")
    rows, results = [], []

    def check(name, value, limit, lower=False):
        ok = math.isfinite(value) and (value >= limit if lower else value <= limit)
        rows.append(dict(check=name, value=value, criterion=(">=" if lower else "<="), threshold=limit,
                         status="PASS" if ok else "FAIL"))
        print(rows[-1], flush=True)
        write_csv(out / "validation.csv", rows)
        if not ok:
            raise RuntimeError("validation failed: " + name)

    def run(name, values, ranks=0, threads=1):
        r = run_case(args.build_dir, args.launcher, out, name, values, ranks, threads, args.timeout)
        results.append(r)
        write_csv(out / "runs.csv", results)
        return r

    base = dict(nx=13, ny=11, nz=9, lx=1, ly=1, lz=0.5, tol=1e-10, steps=20, spinup=10)
    r = run("neumann_discrete", dict(base, mode="poisson-discrete"))
    check("Neumann manufactured discrete solution", r["error_max"], 2e-9)
    check("pressure gauge", abs(r["pressure_mean"]), 1e-10)
    errors = []
    for n in (8, 16, 32):
        r = run("neumann_continuous_"+str(n), dict(base, nx=n, ny=n, nz=n, mode="poisson-continuous"))
        errors.append(r["error_rms"])
    for k in range(2):
        order = math.log(errors[k]/errors[k+1], 2)
        check("pressure spatial order " + str(k+1), order, 1.9, lower=True)
        check("pressure order upper bound " + str(k+1), order, 2.2)

    momentum_errors=[]
    for n in (8,16,32):
        r=run("momentum_"+str(n),dict(base,nx=n,ny=n,nz=n,lx=1,ly=1,lz=1,mode="momentum-test"))
        momentum_errors.append(r["error_rms"])
    for k in range(2):
        order=math.log(momentum_errors[k]/momentum_errors[k+1],2)
        check("momentum spatial order "+str(k+1),order,1.8,lower=True)
        check("momentum order upper bound "+str(k+1),order,2.2)

    def check_projection(r, prefix):
        check(prefix+" divergence reduction", r["div_after_rms"]/r["div_before_rms"], 1.1e-10)
        check(prefix+" residual/projection identity", r["error_max"], 1e-11)
        check(prefix+" zero mean", abs(r["pressure_mean"]), 1e-10)

    r = run("projection_serial", dict(base, mode="projection"))
    check_projection(r, "serial projection")
    r = run("diffusion_cosine", dict(base, force=0, initial="cosine", spinup=0))
    check("3D cosine exact discrete decay", r["pigment_discrete_error"], 2e-12)
    check("zero forcing zero velocity", r["w_max"], 1e-14)

    constant = out / "constant.csv"
    r = run("constant_pigment", dict(base, initial="constant", output=constant))
    check("constant pigment preservation", max(abs(r["c_min"]-1), abs(r["c_max"]-1)), 5e-10)
    check("constant pigment mass", r["mass_abs_error"], 1e-12)

    uniform = out / "uniform.csv"
    r = run("uniform_depth", dict(base, profile="uniform", output=uniform))
    check("uniform forcing vertical velocity", r["w_max"], 1e-11)
    data = read_field(uniform)
    columns = {}
    spread = 0
    for (x,y,z), row in data.items():
        uvp = [row[k] for k in ("u_e", "v_n", "phi")]
        if (x,y) in columns:
            spread = max(spread, max(abs(a-b) for a,b in zip(columns[x,y],uvp)))
        else:
            columns[x,y] = uvp
    check("uniform forcing depth independence of flow", spread, 2e-10)

    serial_path = out / "surface_serial.csv"
    r = run("surface_serial", dict(base, output=serial_path, trace=out / "surface_trace.csv"))
    check("surface forcing produces vertical flow", r["w_max"], 1e-5, lower=True)
    check("surface pigment conservation", r["mass_abs_error"]/abs(r["mass_initial"]), 1e-11)
    check("surface nonnegative pigment", -r["c_min"], 1e-10)
    check("surface bounded pigment", r["c_max"], 1+1e-10)
    check("surface divergence maximum", r["max_divergence"], 1e-8)

    def compare(name, serial_file, parallel_file):
        a,b = read_field(serial_file),read_field(parallel_file)
        if a.keys()!=b.keys():
            raise RuntimeError("missing or mismatched grid cells: " + name)
        for field in ("c","phi","u_e","v_n","w_t"):
            delta = max(abs(a[key][field]-b[key][field]) for key in a)
            check(name+" "+field, delta, 1e-8 if field=="phi" else 2e-10)

    if not args.serial_only:
        layouts = [
            ("mpi1", 1, 1, "blocks", (1,1,1)),
            ("omp2", 1, 2, "blocks", (1,1,1)),
            ("slabs4", 4, 1, "slabs", (1,4,1)),
            ("pencils4", 4, 1, "pencils", (2,2,1)),
            ("blocks4_omp2", 4, 2, "blocks", (2,1,2)),
            ("blocks6", 6, 1, "blocks", (2,1,3)),
        ]
        for name, ranks, threads, decomp, dims in layouts:
            path = out / (name+".csv")
            r = run(name, dict(base, decomp=decomp, px=dims[0], py=dims[1], pz=dims[2], output=path), ranks, threads)
            compare(name, serial_path, path)
        r = run("projection_mpi", dict(base, mode="projection", px=2, py=1, pz=2), 4, 2)
        check_projection(r, "MPI projection")
        r=run("momentum_mpi",dict(base,nx=16,ny=16,nz=16,lx=1,ly=1,lz=1,mode="momentum-test",px=2,py=1,pz=2),4,2)
        check("distributed momentum analytical error agreement",abs(r["error_rms"]-momentum_errors[1]),1e-12)
        # With subdomains only one cell thick, edge ownership and corner exchange
        # cannot hide behind a large interior. Also tests splitting all 3 axes.
        thin = dict(base, nx=3, ny=3, nz=2, steps=10, spinup=5)
        thin_serial,thin_mpi=out/"thin_serial.csv",out/"thin_mpi.csv"
        run("thin_serial", dict(thin, output=thin_serial))
        run("thin_mpi", dict(thin, px=2, py=2, pz=2, output=thin_mpi), 8, 1)
        compare("thin blocks", thin_serial, thin_mpi)
        uniform_mpi = out / "uniform_mpi.csv"
        r = run("uniform_mpi", dict(base, profile="uniform", px=2, py=1, pz=2, output=uniform_mpi), 4, 2)
        compare("uniform depth split", uniform, uniform_mpi)
        check("MPI uniform vertical velocity", r["w_max"], 1e-11)
    (out/"status.json").write_text(json.dumps(dict(status="PASS", tests=len(rows), mpi_tested=not args.serial_only), indent=2)+"\n")
    print("PASS: {} checks; {}".format(len(rows), out/"validation.csv"), flush=True)


if __name__=="__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        sys.exit(1)
