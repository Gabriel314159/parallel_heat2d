#!/usr/bin/env python3
"""Run one performance experiment of experiments/ and record every run.

Expands the file as experiments/README.md ("Run expansion") describes, runs each configuration
warmup_runs + repetitions times, and appends one row per run (failed ones too) to its metrics_csv.

Example:
    python3 scripts/run_benchmarks.py --config experiments/strong_scaling_mpi.json --dry-run
"""

from __future__ import annotations

import argparse
import copy
import csv
import datetime as dt
import hashlib
import itertools
import json
import math
import os
import platform
import re
import shlex
import signal
import shutil
import subprocess
import sys
from pathlib import Path

from heat2d_tools import REPO, deep_merge, load_json, parse_result, problem_args, resolve_base

COLUMNS_FILE = REPO / "results/raw/metrics_columns.csv"
DECOMPOSITIONS = ("rows", "2d")
# Backend -> executable.
EXECUTABLES = {"serial": "heat_serial", "serial_omp": "heat_serial_omp", "parallel": "heat_parallel"}
EXCHANGES = ("blocking", "overlap")
# Columns taken from the configuration, so that error rows have them too.
PROBLEM_COLUMNS = ("nx", "ny", "lx", "ly", "steps", "kappa", "viscosity", "wind", "wind_amplitude",
                   "wind_alpha", "wind_period", "case")
# Placement and binding variables, recorded in session.json.
RECORDED_ENVIRONMENT = ("OMP_PLACES", "OMP_PROC_BIND", "OMP_DISPLAY_AFFINITY", "OMP_NUM_THREADS",
                        "PBS_JOBID", "PBS_QUEUE", "PBS_NODEFILE", "SLURM_JOB_ID", "LOADEDMODULES",
                        "UCX_TLS", "OMPI_MCA_pml", "OMPI_MCA_btl")
PLACEHOLDER = re.compile(r"\{(\w+)\}")


def resolve_path(path: str | Path) -> Path:
    """Paths in experiment files are relative to the repository."""
    p = Path(path).expanduser()
    return p if p.is_absolute() else REPO / p


def metric_columns() -> list[str]:
    with COLUMNS_FILE.open(newline="") as handle:
        return next(csv.reader(handle))


def check_execution(name: str, execution: dict) -> None:
    ranks, threads = execution.get("mpi_ranks"), execution.get("omp_threads")
    if any(type(x) is not int or x < 1 for x in (ranks, threads)):
        raise ValueError(f"{name}: mpi_ranks and omp_threads must be positive integers")
    if execution["backend"] not in EXECUTABLES:
        raise ValueError(f"{name}: backend must be one of {sorted(EXECUTABLES)}, not {execution['backend']!r}")
    if execution["backend"] != "parallel":
        if (ranks, threads) != (1, 1):
            raise ValueError(f"{name}: {EXECUTABLES[execution['backend']]} runs one rank and one thread")
        return
    if execution["decomposition"] not in DECOMPOSITIONS:
        raise ValueError(f"{name}: decomposition must be one of {DECOMPOSITIONS}")
    if execution["exchange"] not in EXCHANGES:
        raise ValueError(f"{name}: exchange must be one of {EXCHANGES}")
    grid = execution.get("process_grid")
    if grid is not None:
        if execution["decomposition"] != "2d":
            raise ValueError(f"{name}: a process_grid needs the 2d decomposition")
        py, px = grid
        if py * px != ranks:
            raise ValueError(f"{name}: process grid {py} x {px} does not match {ranks} ranks")


def configurations(experiment: dict, base: dict) -> list[tuple[str, dict]]:
    """The named configurations of an experiment file, in the order they run."""
    kind = experiment["kind"]
    if kind == "correctness":
        raise ValueError("experiments/correctness.json is run by scripts/check_correctness.py")
    if kind == "strategies":
        combinations = list(itertools.product(experiment["resources"], experiment["decompositions"],
                                              experiment["exchanges"]))
    elif kind in ("strong_scaling", "weak_scaling"):
        combinations = [(r, experiment["decomposition"], experiment["exchange"]) for r in experiment["resources"]]
    else:
        raise ValueError(f"unknown experiment kind: {kind}")

    runs = []
    for resource, decomposition, exchange in combinations:
        if {"decomposition", "exchange"} & set(resource):
            raise ValueError("set decomposition and exchange for the whole file, not per resource")
        execution = {key: value for key, value in resource.items() if key != "problem"}
        config = deep_merge(base, {"execution": {**execution, "decomposition": decomposition, "exchange": exchange},
                                   "problem": resource.get("problem", {})})
        if kind == "weak_scaling":
            derived = {"nx", "ny", "lx", "ly"} & set(resource.get("problem", {}))
            if derived:
                raise ValueError(f"weak scaling derives {sorted(derived)} from the process grid; do not set them")
            py, px = resource["process_grid"]
            problem = config["problem"]
            problem["nx"], problem["ny"] = experiment["local_nx"] * px, experiment["local_ny"] * py
            problem["lx"] = problem["nx"] * experiment["spacing"]["dx"]
            problem["ly"] = problem["ny"] * experiment["spacing"]["dy"]
        e = config["execution"]
        name = f"p{e['mpi_ranks']}_t{e['omp_threads']}_{decomposition}_{exchange}"
        check_execution(name, e)
        if kind == "strategies" and e["mpi_ranks"] * e["omp_threads"] != experiment["total_cores"]:
            raise ValueError(f"{name}: uses {e['mpi_ranks'] * e['omp_threads']} cores, "
                             f"not total_cores = {experiment['total_cores']}")
        problem_args(config["problem"])  # rejects unknown problem keys
        runs.append((name, config))
    names = [name for name, _ in runs]
    if len(set(names)) != len(names):
        raise ValueError("two configurations of this file have the same name")

    if experiment.get("include_serial_reference"):
        if kind == "weak_scaling":
            first = runs[0][1]["execution"]
            if (first["mpi_ranks"], first["omp_threads"]) != (1, 1):
                raise ValueError("the serial reference of weak scaling needs a first entry of one rank and one thread")
            reference = copy.deepcopy(runs[0][1])
        else:
            reference = copy.deepcopy(base)
        for index, (name, backend) in enumerate((("serial_reference", "serial"),
                                                  ("serial_omp_reference", "serial_omp"))):
            config = copy.deepcopy(reference)
            config["execution"] = {"backend": backend, "mpi_ranks": 1, "omp_threads": 1}
            runs.insert(index, (name, config))
    return runs


def placement(execution: dict, args: argparse.Namespace) -> dict[str, int]:
    """Launcher placeholder values defined for this configuration."""
    ranks, threads = execution["mpi_ranks"], execution["omp_threads"]
    values = {"ranks": ranks, "threads": threads}
    if args.cores_per_node:
        if threads > args.cores_per_node:
            raise ValueError(f"{threads} threads of one rank do not fit on {args.cores_per_node} cores of a node")
        nodes = math.ceil(ranks * threads / args.cores_per_node)
        values["nodes"] = nodes
        if ranks % nodes == 0:
            values["ranks_per_node"] = ranks // nodes
            if args.sockets_per_node and ranks % (nodes * args.sockets_per_node) == 0:
                values["ranks_per_socket"] = ranks // (nodes * args.sockets_per_node)
    return values


def fill(token: str, values: dict[str, int], name: str) -> str:
    def value(match: re.Match) -> str:
        key = match.group(1)
        if key not in values:
            raise ValueError(f"{name}: the launcher placeholder {{{key}}} is not defined for "
                             f"{values['ranks']} ranks x {values['threads']} threads (see --help)")
        return str(values[key])
    return PLACEHOLDER.sub(value, token)


def command_for(name: str, config: dict, args: argparse.Namespace) -> list[str]:
    """The command of one configuration, without --output."""
    problem, execution = config["problem"], config["execution"]
    serial = execution["backend"] != "parallel"
    binary = Path(args.build_dir) / EXECUTABLES[execution["backend"]]
    command = [str(binary), *problem_args(problem)]
    if not serial:
        command += ["--decomp", execution["decomposition"], "--exchange", execution["exchange"],
                    "--threads", str(execution["omp_threads"])]
        if execution.get("process_grid"):
            py, px = execution["process_grid"]
            command += ["--py", str(py), "--px", str(px)]
    if serial and not args.launch_serial:
        return command
    values = placement(execution, args)
    return [args.launcher, *(fill(token, values, name) for token in args.mpi_arg),
            "-n", str(execution["mpi_ranks"]), *command]


def head_from_files() -> str | None:
    """HEAD read from .git, for nodes without the git program."""
    git_dir = REPO / ".git"
    try:
        head = (git_dir / "HEAD").read_text().strip()
        if not head.startswith("ref: "):
            return head
        ref = head[5:]
        if (git_dir / ref).is_file():
            return (git_dir / ref).read_text().strip()
        for line in (git_dir / "packed-refs").read_text().splitlines():
            if line.endswith(" " + ref):
                return line.split()[0]
    except OSError:
        pass
    return None


def git_commit() -> str:
    """HEAD, with "-dirty" if tracked files differ from it; "-unchecked" if git is missing and only .git was read."""
    if shutil.which("git") is None:
        head = head_from_files()
        return "unknown" if head is None else head + "-unchecked"
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True)
    if head.returncode != 0:
        return "unknown"
    status = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPO,
                            capture_output=True, text=True)
    return head.stdout.strip() + ("-dirty" if status.stdout.strip() else "")


def environment_record(build_dir: Path, commit: str) -> dict:
    record = {"platform": platform.platform(), "hostname": platform.node(), "python": platform.python_version(),
              "logical_cpus": os.cpu_count(), "git_commit": commit,
              "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "environment": {key: os.environ[key] for key in RECORDED_ENVIRONMENT if key in os.environ}}
    nodefile = os.environ.get("PBS_NODEFILE")
    if nodefile and Path(nodefile).is_file():
        record["pbs_nodes"] = sorted(set(Path(nodefile).read_text().split()))
    if hasattr(os, "sched_getaffinity"):
        record["cpu_affinity"] = sorted(os.sched_getaffinity(0))
    record["binary_sha256"] = {}
    for name in EXECUTABLES.values():
        path = build_dir / name
        if path.is_file():
            record["binary_sha256"][name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return record


def save_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_once(command: list[str], threads: int, timeout: float, log_prefix: Path) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE")
    proc = subprocess.Popen(command, cwd=REPO, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            start_new_session=True)
    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        # Stop the launcher and everything it started.
        timed_out = True
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            stdout, stderr = proc.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            stdout, stderr = proc.communicate()
    log_prefix.with_suffix(".stdout").write_text(stdout)
    log_prefix.with_suffix(".stderr").write_text(stderr)
    if timed_out:
        raise RuntimeError(f"no result after {timeout:g} s (--timeout): {shlex.join(command)}")
    if proc.returncode != 0:
        raise RuntimeError(f"exit code {proc.returncode}: {shlex.join(command)}: {stderr.strip()[-2000:]}")
    return parse_result(stdout)


def check_result(result: dict, config: dict, nodes: int | None, args: argparse.Namespace) -> None:
    """Fail the run unless it ran what was asked and its numbers are sane."""
    problem, execution = config["problem"], config["execution"]
    if result["build_type"] not in ("Release", "RelWithDebInfo"):
        raise RuntimeError(f"benchmarks need a Release or RelWithDebInfo build, not {result['build_type']}")
    asked = {"mpi_ranks": execution["mpi_ranks"], "omp_threads": execution["omp_threads"],
             "nx": problem["nx"], "ny": problem["ny"], "steps": problem["steps"]}
    if execution["backend"] == "parallel":
        asked.update(decomposition=execution["decomposition"], exchange=execution["exchange"])
    if nodes is not None:
        asked["nodes"] = nodes
    wrong = {key: (value, result.get(key)) for key, value in asked.items() if result.get(key) != value}
    if wrong:
        raise RuntimeError("the run differs from the request (asked, used): " +
                           ", ".join(f"{key} {a} vs {u}" for key, (a, u) in wrong.items()))
    for key in ("elapsed_s", "mass_error"):
        value = result.get(key)
        if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise RuntimeError(f"invalid {key}: {value}")
    if result["mass_error"] > 1e-8:
        raise RuntimeError(f"mass not conserved: relative change {result['mass_error']:.2e}")
    if args.require_binding and result["binding"] != "cores":
        raise RuntimeError(f"binding is {result['binding']!r}, not 'cores' (hpc/README.md, \"Binding\")")
    if args.require_gpu and execution["backend"] == "parallel" and result.get("poisson_device") != "gpu":
        raise RuntimeError(f"the Poisson solve ran on {result.get('poisson_device')!r}, not on a GPU")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", "--experiment", dest="config", required=True,
                        help="experiment file (experiments/*.json)")
    parser.add_argument("--repetitions", type=int, help="override the file's repetitions")
    parser.add_argument("--warmup-runs", type=int, help="override the file's warmup_runs")
    parser.add_argument("--run", action="append", metavar="NAME",
                        help="only this configuration (repeatable; --dry-run lists them)")
    parser.add_argument("--dry-run", action="store_true", help="write and print the plan, run nothing")
    parser.add_argument("--keep-going", action="store_true",
                        help="repeat a configuration even after one of its runs failed (default: go to the next)")
    parser.add_argument("--build-dir", default=str(REPO / "build"))
    parser.add_argument("--launcher", default="mpiexec", help="MPI launcher (default mpiexec)")
    parser.add_argument("--mpi-arg", action="append", default=[], metavar="ARG",
                        help="launcher option, repeatable as --mpi-arg=VALUE; may use {ranks} and {threads}")
    parser.add_argument("--launch-serial", action="store_true",
                        help="start the serial references through the launcher too")
    parser.add_argument("--cores-per-node", type=int, help="cores used per node; defines {nodes} and {ranks_per_node}")
    parser.add_argument("--sockets-per-node", type=int, help="with --cores-per-node: defines {ranks_per_socket}")
    parser.add_argument("--allocated-cores", type=int,
                        help="cores of the allocation (default: this process's CPU affinity)")
    parser.add_argument("--require-binding", action="store_true", help="fail runs whose binding is not 'cores'")
    parser.add_argument("--require-gpu", action="store_true",
                        help="fail parallel runs whose Poisson solve did not run on a GPU")
    parser.add_argument("--metrics-csv", help="append to this CSV instead of the file's metrics_csv")
    parser.add_argument("--timeout", type=float, default=3600.0, help="seconds per solver run (default 3600)")
    args = parser.parse_args()
    sys.stdout.reconfigure(line_buffering=True)  # progress in order in batch logs
    args.build_dir = str(Path(args.build_dir).resolve())

    experiment_path = resolve_path(args.config).resolve()
    experiment = load_json(experiment_path)
    if experiment.get("schema_version") != 5:
        raise ValueError(f"{experiment_path.name}: expected schema_version 5")
    base = resolve_base(experiment_path)
    runs = configurations(experiment, base)
    if args.run:
        unknown = set(args.run) - {name for name, _ in runs}
        if unknown:
            raise ValueError(f"no such configuration in {experiment_path.name}: {sorted(unknown)}")
        runs = [(name, config) for name, config in runs if name in args.run]
    repetitions = experiment["repetitions"] if args.repetitions is None else args.repetitions
    warmups = experiment["warmup_runs"] if args.warmup_runs is None else args.warmup_runs
    if repetitions < 1 or warmups < 0:
        raise ValueError("repetitions must be positive and warmup runs non-negative")
    if args.sockets_per_node and not args.cores_per_node:
        raise ValueError("--sockets-per-node needs --cores-per-node")

    # All commands first, so that a bad placeholder stops the plan here.
    commands = {name: command_for(name, config, args) for name, config in runs}
    session = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = resolve_path(base["output"]["directory"]) / session
    out.mkdir(parents=True, exist_ok=False)
    csv_path = resolve_path(args.metrics_csv or experiment["metrics_csv"])
    commit = git_commit()
    metadata = environment_record(Path(args.build_dir), commit)
    metadata.update(config_file=str(experiment_path), kind=experiment["kind"], session_id=session,
                    metrics_csv=str(csv_path), options=vars(args), repetitions=repetitions, warmup_runs=warmups,
                    runs=[{"name": name, "command": commands[name], **config} for name, config in runs])
    save_json(out / "session.json", metadata)
    for name, _ in runs:
        print(f"{name}: {shlex.join(commands[name])}")
    if args.dry_run:
        print(f"{len(runs)} configurations, {warmups} warmup(s) + {repetitions} repetition(s) each; "
              f"nothing run. Plan: {out / 'session.json'}")
        return 0

    available = args.allocated_cores or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity")
                                         else os.cpu_count() or 1)
    for name, config in runs:
        e = config["execution"]
        if e["mpi_ranks"] * e["omp_threads"] > available:
            raise ValueError(f"{name} needs {e['mpi_ranks'] * e['omp_threads']} cores; {available} are available "
                             "(--allocated-cores)")
    needed = {EXECUTABLES[config["execution"]["backend"]] for _, config in runs}
    for binary in sorted(needed):
        if not (Path(args.build_dir) / binary).is_file():
            raise FileNotFoundError(f"{Path(args.build_dir) / binary} not found; build first (Release)")
    if any(config["execution"]["backend"] == "parallel" or args.launch_serial for _, config in runs):
        if shutil.which(args.launcher) is None:
            raise FileNotFoundError(f"MPI launcher not found: {args.launcher}")

    columns = metric_columns()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    existed = csv_path.is_file() and csv_path.stat().st_size > 0
    if existed:
        with csv_path.open(newline="") as handle:
            if next(csv.reader(handle)) != columns:
                raise ValueError(f"{csv_path} has other columns than {COLUMNS_FILE.name}; choose a new metrics_csv")

    failures = 0
    with csv_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore", restval="")
        if not existed:
            writer.writeheader()
        for name, config in runs:
            problem, execution = config["problem"], config["execution"]
            ranks, threads = execution["mpi_ranks"], execution["omp_threads"]
            for index in range(warmups + repetitions):
                phase = "warmup" if index < warmups else "measured"
                repetition = index if phase == "warmup" else index - warmups
                row = {"experiment": experiment_path.stem, "session_id": session, "run_id": name,
                       "phase": phase, "repetition": repetition, "status": "ok", "backend": execution["backend"],
                       "mpi_ranks": ranks, "omp_threads": threads, "total_cores": ranks * threads,
                       "hostname": platform.node(), "git_commit": commit,
                       "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
                row.update({key: problem.get(key) for key in PROBLEM_COLUMNS})
                prefix = out / f"{name}_{phase}_{repetition}"
                command = list(commands[name])
                if config["output"].get("write_field"):
                    command += ["--output", str(prefix.with_suffix(".csv"))]
                nodes = placement(execution, args).get("nodes") if execution["backend"] == "parallel" else None
                try:
                    result = run_once(command, threads, args.timeout, prefix)
                    check_result(result, config, nodes, args)
                    row.update({key: value for key, value in result.items() if key in columns})
                    print(f"{name}: {phase} {repetition}, {result['elapsed_s']:.4g} s, "
                          f"{result['poisson_iterations_mean']:.0f} CG iterations per step, binding {result['binding']}")
                except Exception as error:  # any failure becomes an error row
                    row.update(status="error", error=" | ".join(str(error).replace("\0", "").splitlines()))
                    failures += 1
                    print(f"{name}: {phase} {repetition} FAILED: {error}", file=sys.stderr)
                writer.writerow(row)
                handle.flush()
                if row["status"] == "error" and not args.keep_going:
                    break  # the rest of this configuration; the next one still runs

    print(f"Raw measurements: {csv_path}")
    print(f"Resolved plan, environment and logs: {out}")
    if failures:
        print(f"{failures} run(s) failed; their rows have status 'error'", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, FileNotFoundError, OSError) as error:
        print(f"run_benchmarks: {error}", file=sys.stderr)
        sys.exit(2)
