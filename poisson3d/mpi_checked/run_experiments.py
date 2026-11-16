#!/usr/bin/env python3
"""Validate the 3D Poisson solver or run its Vanda experiments."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

LOG_DIRECTORY = None
LAUNCH_NUMBER = 0


@dataclass(frozen=True)
class Case:
    experiment: str
    grid: tuple[int, int, int]
    lengths: tuple[float, float, float] = (1.0, 1.0, 1.0)
    ranks: int = 1
    threads: int = 1
    process_grid: tuple[int, int] | None = None  # px, py
    backend: str = "mpi"
    exchange: str = "blocking"
    nodes_override: int | None = None

    @property
    def nodes(self):
        if self.nodes_override is not None:
            return self.nodes_override
        if self.experiment == "openmp" or self.backend in ("cpu-omp", "gpu"):
            return 1
        return (self.ranks * self.threads + 63) // 64


EXECUTABLES = {"serial": "poisson3d_serial", "serial-omp": "poisson3d_serial_omp", "mpi": "poisson3d_mpi",
               "cpu-omp": "poisson3d_omp", "gpu": "poisson3d_gpu"}


def fixed(experiment, n, ranks=1, threads=1, process_grid=None, backend="mpi"):
    return Case(experiment, (n, n, n), ranks=ranks, threads=threads,
                process_grid=process_grid, backend=backend)


def weak(ranks):
    p = math.isqrt(ranks)
    n = p * 32
    return Case("weak", (n, n, 32), ((n + 1) / 33, (n + 1) / 33, 1.0),
                ranks=ranks, process_grid=(p, p))


def plan(suite, one_gpu_only=False, include_controls=False):
    if suite == "one-node":
        return (
            [fixed("serial-size", n, backend="serial") for n in (32, 48, 64, 96, 128)]
            + [fixed("strong", 128, backend=b) for b in ("serial", "serial-omp")]
            + [fixed("strong", 128, p) for p in (1, 2, 4, 8, 16, 32, 64)]
            + [fixed("openmp", 128, backend=b) for b in ("serial", "serial-omp")]
            + [fixed("openmp", 128, threads=t) for t in (1, 2, 4, 8, 16, 32, 36, 64, 72)]
            + ([fixed("hybrid-node", 128, 64 // t, t) for t in (1, 2, 4, 8, 16, 32)] if include_controls else [])
        )
    if suite == "four-nodes":
        # Four rows per rank for row decomposition at 256 ranks, so some
        # interior rows can be computed while messages are in flight.
        communication = [
            Case("communication", (256, 1024, 32), (1.0, 1025 / 257, 33 / 257),
                 ranks=256, process_grid=shape, exchange=exchange)
            for shape in ((1, 256), (16, 16)) for exchange in ("blocking", "overlap")
        ]
        spread = [
            Case("node-spread", grid, lengths, ranks=64, process_grid=(8, 8), nodes_override=nodes)
            for grid, lengths in [((128, 128, 128), (1.0, 1.0, 1.0)), (weak(64).grid, weak(64).lengths)]
            for nodes in (1, 2, 4)
        ]
        return (
            [fixed("strong", 128, p) for p in (128, 256)]
            + [fixed("weak", 32, backend=b) for b in ("serial", "serial-omp")]
            + [weak(p) for p in (1, 4, 9, 16, 36, 64, 144, 256)]
            + [fixed("hybrid", 128, 256 // t, t) for t in (1, 2, 4, 8, 16, 32, 64)]
            + communication + (spread if include_controls else [])
        )
    if suite == "gpu":
        return [fixed("gpu-size", n, ranks=p, threads=36, backend=backend)
                for n in (64, 128, 192, 256)
                for p in ((1,) if one_gpu_only else (1, 2))
                for backend in ("cpu-omp", "gpu")]
    raise ValueError(f"unexpected suite {suite}")


def placement(case):
    if case.experiment == "openmp" and case.backend == "mpi":
        return "{0}:36:2,{1}:36:2"
    if case.experiment == "hybrid" and case.threads == 64:
        return "{0}:32:2,{1}:32:2"
    return "cores"


def launcher(case, *, verify=False):
    args = shlex.split(os.environ.get("POISSON3D_MPIRUN", "mpirun"))
    if os.environ.get("EBROOTPRRTE"):
        args += ["--prefix", os.environ["EBROOTPRRTE"]]
    if os.environ.get("LD_LIBRARY_PATH"):
        args += ["-x", "LD_LIBRARY_PATH"]
    args += ["-n", str(case.ranks), "-x", "OMP_NUM_THREADS",
             "-x", "OMP_PLACES", "-x", "OMP_PROC_BIND"]
    if case.backend == "gpu":
        args += ["-x", "OMP_TARGET_OFFLOAD"]
        if "CUDA_VISIBLE_DEVICES" in os.environ:
            args += ["-x", "CUDA_VISIBLE_DEVICES"]
    args += ["--report-bindings"]
    if case.experiment == "preflight":
        return args + ["--map-by", "ppr:1:node", "--bind-to", "core"]
    if verify:
        return args + ["--oversubscribe", "--bind-to", "none"]
    if placement(case) != "cores":
        return args + ["--map-by", "ppr:1:node", "--bind-to", "none"]
    if case.threads == 1 and case.ranks <= 32:
        return args + ["--map-by", "core", "--bind-to", "core"]
    if case.experiment == "node-spread":
        socket_cores = case.ranks // (2 * case.nodes)
    elif case.experiment == "weak" and case.ranks % (2 * case.nodes) == 0:
        # Match the 2D follow-up: 36 -> 18/socket, 144 -> 24/socket.
        socket_cores = case.ranks // (2 * case.nodes)
    else:
        socket_cores = 36 if case.backend in ("cpu-omp", "gpu") else 32
    return args + ["--map-by", f"ppr:{socket_cores // case.threads}:socket:PE={case.threads}",
                   "--bind-to", "core"]


def command(case, build, *, timed=False, rhs="benchmark", tol="1e-8", solution=None):
    executable = EXECUTABLES[case.backend]
    args = [str(build / executable)]
    for axis, n, length in zip("xyz", case.grid, case.lengths):
        args += [f"--n{axis}", str(n), f"--l{axis}", str(length)]
    if case.process_grid:
        px, py = case.process_grid
        args += ["--px", str(px), "--py", str(py)]
    args += ["--rhs", rhs, "--tol", tol, "--exchange", case.exchange]
    if solution is not None:
        args += ["--solution", str(solution)]
    return (launcher(case, verify=not timed) if timed or case.backend in ("mpi", "cpu-omp", "gpu") else []) + args


def run(case, cmd, timeout, *, timed=False):
    env = dict(os.environ, OMP_NUM_THREADS=str(case.threads), OMP_PLACES=placement(case) if timed else "cores",
               OMP_PROC_BIND="close" if timed else "false")
    if case.backend == "gpu":
        env["OMP_TARGET_OFFLOAD"] = "MANDATORY"
    proc = execute(cmd, env, timeout)
    lines = [line for line in proc.stdout.splitlines() if line.startswith('{"nx":')]
    if len(lines) != 1:
        raise RuntimeError(f"expected a single result from {shlex.join(cmd)}\n{proc.stdout}\n{proc.stderr}")
    out = json.loads(lines[0])
    if (out["ranks"], out["threads"]) != (case.ranks, case.threads):
        raise RuntimeError("observed MPI/OpenMP layout differs from the requested layout")
    expected_device = "gpu" if case.backend == "gpu" else "cpu"
    if out.get("poisson_device", "cpu") != expected_device:
        raise RuntimeError(f"expected solver device {expected_device}, got {out}")
    if out["gpu_count"] != (case.ranks if case.backend == "gpu" else 0):
        raise RuntimeError(f"observed GPU count does not match the requested layout: {out}")
    if timed and (out["binding"] != "cores" or out["nodes"] != case.nodes):
        raise RuntimeError(f"unexpected placement: {out}. Check MPI mapping and the PBS allocation.")
    return out


def execute(cmd, env, timeout):
    global LAUNCH_NUMBER
    LAUNCH_NUMBER += 1
    prefix = LOG_DIRECTORY / f"launch-{LAUNCH_NUMBER:04d}" if LOG_DIRECTORY else None
    print(f"START launch {LAUNCH_NUMBER}: {shlex.join(cmd)}", flush=True)
    if prefix:
        prefix.with_suffix(".json").write_text(json.dumps({
            "command": cmd, "timeout": timeout,
            "environment": {k: env.get(k) for k in
                            ("OMP_NUM_THREADS", "OMP_PLACES", "OMP_PROC_BIND", "EBROOTPRRTE", "LD_LIBRARY_PATH")},
        }, indent=2) + "\n")
    try:
        proc = subprocess.run(cmd, env=env, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        for kind in ("stdout", "stderr"):
            value = getattr(exc, kind) or ""
            if isinstance(value, bytes):
                value = value.decode(errors="replace")
            if prefix:
                prefix.with_suffix(f".{kind}").write_text(value)
        raise RuntimeError(f"MPI launch timed out after {timeout}s: {shlex.join(cmd)}. See {prefix}") from exc
    if prefix:
        prefix.with_suffix(".stdout").write_text(proc.stdout)
        prefix.with_suffix(".stderr").write_text(proc.stderr)
    if proc.returncode:
        raise RuntimeError(f"exit {proc.returncode}: {shlex.join(cmd)}\n{proc.stdout}\n{proc.stderr}")
    return proc


def preflight(build, hosts, timeout):
    """Check remote launch, topology, then actual MPI traffic before the suite."""
    require(len(hosts) == 4, "four-node preset requires exactly four allocated hosts")
    case = Case("preflight", (23, 29, 19), (1.5, 0.9, 1.2),
                ranks=4, process_grid=(2, 2), nodes_override=4)
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_PLACES="cores", OMP_PROC_BIND="close")
    # Probe without rank binding so it can inspect the entire allocated CPU set.
    cmd = launcher(case)
    cmd[-1] = "none"
    cmd += ["python3", str(Path(__file__).with_name("node_probe.py").resolve())]
    proc = execute(cmd, env, min(timeout, 120))
    records = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith('{"host":')]
    require(len(records) == 4 and len({r["host"] for r in records}) == 4,
            f"expected one probe per node, got {proc.stdout}")
    require({r["host"].split('.')[0].lower() for r in records} ==
            {h.split('.')[0].lower() for h in hosts}, "probe hosts differ from PBS allocation")
    require(all(r["vanda_layout_ok"] for r in records),
            f"nodes do not match the 72-core, two-socket layout assumed by this benchmark: {records}")
    (LOG_DIRECTORY.parent / "preflight-nodes.json").write_text(json.dumps(records, indent=2) + "\n")
    out = run(case, command(case, build, timed=True, rhs="discrete", tol="1e-10"),
              min(timeout, 120), timed=True)
    require(out["max_error"] < 1e-9, f"cross-node numerical preflight failed: {out}")
    (LOG_DIRECTORY.parent / "preflight-solver.json").write_text(json.dumps(out, indent=2) + "\n")
    print("Four-node launch, topology, binding and numerical preflight passed.", flush=True)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def verify(build, timeout, serial_only):
    continuous = []
    for n in (15, 31, 63):
        case = fixed("continuous", n, backend="serial")
        out = run(case, command(case, build, rhs="continuous", tol="1e-10"), timeout)
        continuous.append(out)
        print(f"continuous reference {n}^3: relative L2 error {out['relative_error']:.3e}")
    errors = [r["relative_error"] for r in continuous]
    rates = [math.log2(a / b) for a, b in zip(errors, errors[1:])]
    require(min(rates) > 1.8, f"expected second-order accuracy, got {rates}")
    print("observed convergence orders:", ", ".join(f"{rate:.2f}" for rate in rates))

    # Uneven blocks in x and y on an anisotropic grid exercise both halos.
    grid, lengths = (23, 29, 19), (1.5, 0.9, 1.2)
    with TemporaryDirectory() as tmp:
        ref, actual = Path(tmp) / "serial.txt", Path(tmp) / "parallel.txt"
        case = Case("discrete", grid, lengths, backend="serial")
        out = run(case, command(case, build, rhs="discrete", tol="1e-10", solution=ref), timeout)
        require(out["max_error"] < 1e-9, f"discrete reference check failed: {out}")
        print(f"discrete reference on {grid}: max error {out['max_error']:.3e}")
        validation = dict(continuous=continuous, orders=rates, discrete_serial=out, parallel=[])
        if serial_only:
            return validation
        expected = [float(value) for value in ref.read_text().splitlines()]
        for p, t, shape, exchange in [(1, 2, (1, 1), "blocking"), (3, 1, (3, 1), "blocking"),
                                      (4, 1, (2, 2), "blocking"), (6, 2, (2, 3), "overlap"),
                                      (4, 1, (1, 4), "overlap")]:
            case = Case("discrete", grid, lengths, p, t, shape, exchange=exchange)
            cmd = command(case, build, rhs="discrete", tol="1e-10", solution=actual)
            out = run(case, cmd, timeout)
            require(out["max_error"] < 1e-9, f"discrete reference check failed: {out}")
            field = [float(value) for value in actual.read_text().splitlines()]
            require(len(field) == len(expected), "serial and MPI field sizes differ")
            difference = max(abs(a - b) for a, b in zip(field, expected))
            require(difference < 1e-9, f"serial and MPI fields differ by {difference:.3e}")
            validation["parallel"].append(dict(out, max_serial_parallel_difference=difference))
            print(f"{p} ranks x {t} threads, grid {shape}, {exchange}: max field difference {difference:.3e}")
        return validation


def verify_gpu(build, timeout, one_gpu_only=False):
    grid, lengths = (23, 29, 19), (1.5, 0.9, 1.2)
    validation = []
    with TemporaryDirectory() as tmp:
        shapes = [(1, 1)] if one_gpu_only else [(1, 1), (1, 2), (2, 1)]
        for shape in shapes:
            for rhs in ("discrete", "benchmark"):
                records, fields = [], []
                for backend in ("cpu-omp", "gpu"):
                    case = Case("gpu-validation", grid, lengths, ranks=math.prod(shape),
                                process_grid=shape, backend=backend)
                    path = Path(tmp) / f"{backend}.txt"
                    cmd = command(case, build, rhs=rhs, tol="1e-10", solution=path)
                    out = run(case, cmd, timeout)
                    if rhs == "discrete":
                        require(out["max_error"] < 1e-9, f"discrete reference check failed: {out}")
                    records.append(out)
                    fields.append([float(value) for value in path.read_text().splitlines()])
                require(len(fields[0]) == len(fields[1]), "CPU and GPU fields differ in size")
                difference = max(abs(a - b) for a, b in zip(*fields))
                require(difference < 1e-9, f"CPU and GPU fields differ by {difference:.3e}")
                validation.append(dict(rhs=rhs, cpu=records[0], gpu=records[1], max_difference=difference))
                print(f"GPU validation {shape} ({rhs}): max CPU/GPU difference {difference:.3e}")
    return validation


def benchmark(args, cases):
    global LOG_DIRECTORY
    args.output.mkdir(parents=True, exist_ok=True)
    raw = args.output / f"{args.suite}.csv"
    metadata = args.output / f"{args.suite}.json"
    if raw.exists() or metadata.exists():
        raise FileExistsError(f"choose a fresh output folder, {raw.name} or {metadata.name} already exists")
    LOG_DIRECTORY = args.output / "launch-logs"
    LOG_DIRECTORY.mkdir(exist_ok=False)
    nodefile = os.environ.get("PBS_NODEFILE")
    slots = Path(nodefile).read_text().splitlines() if nodefile else []
    hosts = sorted(set(slots))
    if len(hosts) < max(case.nodes for case in cases):
        raise RuntimeError("run this preset inside its Vanda PBS allocation")
    (args.output / "allocation.json").write_text(json.dumps({
        "pbs_job_id": os.environ.get("PBS_JOBID"), "hosts": hosts, "slots": slots,
        "modules": os.environ.get("LOADEDMODULES"),
    }, indent=2) + "\n")
    if args.suite == "four-nodes":
        preflight(args.build_dir, hosts, args.timeout)
    if args.suite == "gpu":
        validation = verify_gpu(args.build_dir, args.timeout, args.one_gpu_only)
        gpu_inventory = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=name,driver_version,memory.total,pci.bus_id",
             "--format=csv"], text=True).strip()
    else:
        validation = verify(args.build_dir, args.timeout, serial_only=False)
        gpu_inventory = None
    binaries = (("poisson3d_omp", "poisson3d_gpu") if args.suite == "gpu" else
                ("poisson3d_serial", "poisson3d_serial_omp", "poisson3d_mpi"))
    meta = {
        "pbs_job_id": os.environ.get("PBS_JOBID"), "allocated_hosts": hosts, "allocated_slots": len(slots),
        "platform": platform.platform(), "modules": os.environ.get("LOADEDMODULES"),
        "warmup_runs": args.warmups, "repetitions": args.repetitions,
        "binary_sha256": {
            name: hashlib.sha256((args.build_dir / name).read_bytes()).hexdigest()
            for name in binaries
        },
        "gpu_inventory": gpu_inventory,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "commands": [command(case, args.build_dir, timed=True) for case in cases],
        "host_placements": [placement(case) for case in cases],
        "validation": validation,
    }
    metadata.write_text(json.dumps(meta, indent=2) + "\n")
    with raw.open("x", newline="") as handle:
        writer = None
        for case in cases:
            cmd = command(case, args.build_dir, timed=True)
            print(case.experiment, shlex.join(cmd), flush=True)
            for phase, count in (("warmup", args.warmups), ("measured", args.repetitions)):
                for rep in range(1, count + 1):
                    metrics = run(case, cmd, args.timeout, timed=True)
                    row = dict(experiment=case.experiment, phase=phase, repetition=rep,
                               backend=case.backend, **metrics)
                    row["solve_s_per_iteration"] = row["solve_s"] / row["iterations"]
                    if writer is None:
                        writer = csv.DictWriter(handle, fieldnames=list(row))
                        writer.writeheader()
                    writer.writerow(row)
                    handle.flush()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", choices=("verify", "verify-gpu", "one-node", "four-nodes", "gpu"))
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=900)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--serial-only", action="store_true", help="run only the serial validation cases")
    parser.add_argument("--one-gpu-only", action="store_true", help="omit two-GPU configurations")
    parser.add_argument("--include-controls", action="store_true", help="add node-placement controls to a CPU suite")
    args = parser.parse_args()
    args.build_dir = args.build_dir.resolve()
    if args.warmups < 0 or args.repetitions < 1 or args.timeout <= 0:
        parser.error("invalid repetition count or timeout")
    if args.one_gpu_only and args.suite not in ("gpu", "verify-gpu"):
        parser.error("--one-gpu-only is for gpu and verify-gpu")
    if args.include_controls and args.suite not in ("one-node", "four-nodes"):
        parser.error("--include-controls is for one-node and four-nodes")
    if args.suite in ("verify", "verify-gpu"):
        if args.dry_run:
            parser.error("--dry-run is for the benchmark presets")
        if args.suite == "verify-gpu":
            if args.serial_only:
                parser.error("--serial-only is for verify")
            verify_gpu(args.build_dir, args.timeout, args.one_gpu_only)
        else:
            verify(args.build_dir, args.timeout, args.serial_only)
        return
    if args.serial_only:
        parser.error("--serial-only is for verify")
    cases = plan(args.suite, args.one_gpu_only, args.include_controls)
    if args.dry_run:
        for case in cases:
            env = f"OMP_NUM_THREADS={case.threads} OMP_PLACES={shlex.quote(placement(case))} OMP_PROC_BIND=close"
            if case.backend == "gpu":
                env += " OMP_TARGET_OFFLOAD=MANDATORY"
            print(f"{env} {shlex.join(command(case, args.build_dir, timed=True))}")
    else:
        benchmark(args, cases)


if __name__ == "__main__":
    main()
