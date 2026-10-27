#!/usr/bin/env python3
"""Run the correctness cases of experiments/correctness.json and check them.

For every run, the serial solver is run on the resolved problem and checked
(docs/experiment_plan.md, "Correctness"):
  - mass conservation: relative mass change <= mass_conservation_rtol,
  - constant case: the field stays at the background value,
  - cosine_decay: agreement with the exact discrete solution
    (cosine_decay_discrete_rtol); the error against the continuous solution is
    reported separately, as discretization error.

Each run also names a parallel configuration. Once heat_parallel is built,
its final field is compared cell by cell with the serial one
(serial_parallel_atol / serial_parallel_rtol). Until then those comparisons
are reported as skipped; --serial-only skips them silently.

Exit status: 0 if every check passed, 1 otherwise.

Example:
    python3 scripts/check_correctness.py --serial-only
"""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

from heat2d_tools import REPO, deep_merge, executable, load_json, problem_args, read_field_csv, run_solver


def resolve_runs(experiment_path: Path) -> list[dict]:
    """Defaults, then the experiment-level overrides, then each run's overrides."""
    experiment = load_json(experiment_path)
    defaults = load_json((experiment_path.parent / experiment["defaults"]).resolve())
    shared = {key: experiment[key] for key in ("problem", "execution", "output") if key in experiment}
    base = deep_merge(defaults, shared)
    runs = []
    for run in experiment["runs"]:
        resolved = deep_merge(base, {key: run[key] for key in ("problem", "execution", "output") if key in run})
        resolved["name"] = run["name"]
        runs.append(resolved)
    return runs


def parallel_args(execution: dict) -> list[str]:
    args = ["--decomp", execution["decomposition"], "--exchange", execution["exchange"]]
    if execution.get("process_grid"):
        py, px = execution["process_grid"]
        args += ["--py", str(py), "--px", str(px)]
    return args


def max_field_difference(a: list[list[float]], b: list[list[float]], atol: float, rtol: float) -> tuple[float, bool]:
    if len(a) != len(b) or any(len(ra) != len(rb) for ra, rb in zip(a, b)):
        return float("inf"), False
    worst = 0.0
    within = True
    for row_a, row_b in zip(a, b):
        for x, y in zip(row_a, row_b):
            difference = abs(x - y)
            worst = max(worst, difference)
            within = within and difference <= atol + rtol * abs(y)
    return worst, within


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--experiment", default=str(REPO / "experiments/correctness.json"))
    parser.add_argument("--build-dir", default=str(REPO / "build"))
    parser.add_argument("--output-dir", help="where fields are written (default: the experiment's output.directory)")
    parser.add_argument("--serial-only", action="store_true", help="skip the serial/parallel comparisons")
    parser.add_argument("--launcher", default="mpiexec -n {ranks}", help="MPI launch command; {ranks} is replaced")
    parser.add_argument("--run", action="append", help="only this run name (may be repeated)")
    args = parser.parse_args()

    experiment_path = Path(args.experiment).resolve()
    checks = load_json(experiment_path)["checks"]
    runs = resolve_runs(experiment_path)
    if args.run:
        runs = [run for run in runs if run["name"] in args.run]

    serial = executable(args.build_dir, "heat_serial")
    parallel_path = Path(args.build_dir) / "heat_parallel"
    parallel_built = parallel_path.is_file()

    report = []
    failures = 0
    for run in runs:
        name = run["name"]
        problem = run["problem"]
        out_dir = Path(args.output_dir or REPO / run["output"]["directory"])
        out_dir.mkdir(parents=True, exist_ok=True)
        results = []  # (check, passed or None for skipped/info, detail)

        serial_csv = out_dir / f"{name}_serial.csv"
        try:
            summary = run_solver([str(serial), *problem_args(problem), "--output", str(serial_csv)])
        except (RuntimeError, ValueError) as error:
            results.append(("serial run", False, str(error)))
            summary = None

        if summary is not None:
            tol = checks["mass_conservation_rtol"]
            results.append(("mass conservation", summary["mass_error"] <= tol,
                            f"relative change {summary['mass_error']:.2e} (tolerance {tol:.0e})"))

            if problem["case"] == "constant" and checks.get("constant_field_invariance"):
                c0 = problem["background"]
                field = read_field_csv(serial_csv)
                deviation = max(abs(value - c0) for row in field for value in row)
                tol = 1e-12 * max(1.0, abs(c0))
                results.append(("constant field", deviation <= tol, f"max |c - {c0:g}| = {deviation:.2e}"))

            if problem["case"] == "cosine_decay":
                tol = checks["cosine_decay_discrete_rtol"]
                error = summary["cosine_discrete_error"]
                results.append(("cosine vs exact discrete", error <= tol, f"{error:.2e} (tolerance {tol:.0e})"))
                results.append(("cosine vs continuous", None,
                                f"{summary['cosine_continuous_error']:.2e} (discretization error, reported only)"))

            execution = run["execution"]
            if execution.get("backend") == "parallel" and not args.serial_only:
                if not parallel_built:
                    results.append(("serial vs parallel", None, "skipped: heat_parallel is not built yet"))
                else:
                    ranks = int(execution["mpi_ranks"])
                    threads = int(execution["omp_threads"])
                    parallel_csv = out_dir / f"{name}_parallel.csv"
                    command = [*shlex.split(args.launcher.format(ranks=ranks)), str(parallel_path),
                               *problem_args(problem), *parallel_args(execution), "--output", str(parallel_csv)]
                    try:
                        run_solver(command, threads=threads)
                        worst, within = max_field_difference(
                            read_field_csv(parallel_csv), read_field_csv(serial_csv),
                            checks["serial_parallel_atol"], checks["serial_parallel_rtol"])
                        results.append(("serial vs parallel", within,
                                        f"max |difference| = {worst:.2e} ({ranks} ranks x {threads} threads)"))
                    except (RuntimeError, ValueError) as error:
                        results.append(("parallel run", False, str(error)))

        print(f"{name}  ({problem['nx']} x {problem['ny']}, {problem['steps']} steps, case {problem['case']})")
        for check, passed, detail in results:
            mark = "ok  " if passed else ("--  " if passed is None else "FAIL")
            print(f"  {mark}  {check}: {detail}")
            failures += passed is False
        report.append({"run": name, "checks": [
            {"check": check, "passed": passed, "detail": detail} for check, passed, detail in results]})

    report_dir = Path(args.output_dir or REPO / "outputs/correctness")
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"\n{'PASS' if failures == 0 else 'FAIL'}: {failures} failed check(s); report in {report_dir / 'report.json'}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
