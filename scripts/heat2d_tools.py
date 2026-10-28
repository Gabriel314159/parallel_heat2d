"""Shared helpers for the experiment scripts.

- Load and merge the JSON configuration (config/defaults.json, experiments/*.json).
- Turn a resolved configuration into the solvers' command-line flags.
- Run a solver and parse the RESULT line it prints at the end.

Standard library only, so the correctness checks run wherever Python 3 does.
"""

from __future__ import annotations

import copy
import csv
import json
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# JSON problem key -> solver flag (config/README.md). blob_centre is a pair.
PROBLEM_FLAGS = {
    "nx": "--nx",
    "ny": "--ny",
    "lx": "--lx",
    "ly": "--ly",
    "steps": "--steps",
    "kappa": "--kappa",
    "cfl": "--cfl",
    "dt": "--dt",
    "max_speed": "--max-speed",
    "check_interval": "--check-interval",
    "viscosity": "--viscosity",
    "wind": "--wind",
    "wind_amplitude": "--wind-amplitude",
    "wind_alpha": "--wind-alpha",
    "wind_period": "--wind-period",
    "spinup": "--spinup",
    "case": "--case",
    "background": "--background",
    "amplitude": "--amplitude",
    "blob_width": "--blob-width",
}


def load_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def deep_merge(base: dict, override: dict) -> dict:
    """Merge override into a copy of base: nested objects merge, scalars and lists replace."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def format_value(value) -> str:
    """Numbers in a form the C++ parser reads back exactly (repr round-trips floats)."""
    if isinstance(value, bool):
        raise ValueError("boolean values have no solver flag")
    if isinstance(value, float):
        return repr(value)
    return str(value)


def problem_args(problem: dict) -> list[str]:
    """Solver flags for a resolved "problem" object."""
    unknown = set(problem) - set(PROBLEM_FLAGS) - {"blob_centre"}
    if unknown:
        raise ValueError(f"unknown problem keys: {sorted(unknown)}")
    args: list[str] = []
    for key, flag in PROBLEM_FLAGS.items():
        value = problem.get(key)
        if value is not None:
            args += [flag, format_value(value)]
    if problem.get("blob_centre") is not None:
        x, y = problem["blob_centre"]
        args += ["--blob-x", format_value(x), "--blob-y", format_value(y)]
    return args


def executable(build_dir: str | Path, name: str) -> Path:
    path = Path(build_dir) / name
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; build first: cmake -S . -B build && cmake --build build")
    return path


def parse_result(stdout: str) -> dict:
    """The JSON object on the solver's RESULT line."""
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT "):
            return json.loads(line[len("RESULT "):])
    raise ValueError("solver output has no RESULT line")


def run_solver(command: list[str], threads: int = 1, timeout: float | None = None) -> dict:
    """Run one solver command and return its RESULT. Raises RuntimeError on failure."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    completed = subprocess.run(command, cwd=REPO, env=env, capture_output=True, text=True, timeout=timeout)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {completed.returncode}: {' '.join(command)}\n{completed.stderr.strip()}"
        )
    return parse_result(completed.stdout)


def read_field_csv(path: str | Path) -> list[list[float]]:
    with open(path, newline="") as handle:
        return [[float(value) for value in row] for row in csv.reader(handle) if row]
