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
    "poisson_tol": "--poisson-tol",
    "poisson_max_iter": "--poisson-max-iter",
    # The 3D pigment model (heat3d_serial, heat3d_parallel; docs/design_3d.md, "Model A").
    "nz": "--nz",
    "lz": "--lz",
    "blob_z": "--blob-z",  # a depth, or "column"
}
# JSON problem key -> solver flag given without a value when the key is true.
SWITCH_FLAGS = {
    "freeze_flow": "--freeze-flow",
}

# The two models, by the "model" key of an experiment file (default "2d").
MODELS = {
    "2d": {"serial": "heat_serial", "serial_omp": "heat_serial_omp", "parallel": "heat_parallel",
           "decompositions": ("rows", "2d")},
    "3d": {"serial": "heat3d_serial", "parallel": "heat3d_parallel", "decompositions": ("slabs", "pencils")},
}


def model_of(experiment: dict) -> str:
    """The model an experiment file is for: "2d" (default) or "3d"."""
    model = experiment.get("model", "2d")
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; use one of {sorted(MODELS)}")
    return model


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


# The objects of an experiment file that are merged on top of its defaults.
SHARED_KEYS = ("problem", "execution", "output")


def resolve_base(path: str | Path, _depth: int = 0) -> dict:
    """The configuration that a JSON file starts from.

    A file may name another with "defaults" (a path relative to the file). That
    file is resolved first, in the same way, and the problem, execution and
    output objects of this file are merged on top. An experiment can therefore
    inherit experiments/benchmark_problem.json, which inherits config/defaults.json.
    """
    path = Path(path).resolve()
    data = load_json(path)
    if "defaults" not in data:
        return data
    if _depth >= 4:
        raise ValueError(f"{path}: the chain of 'defaults' files is too long (a cycle?)")
    base = resolve_base(path.parent / data["defaults"], _depth + 1)
    return deep_merge(base, {key: data[key] for key in SHARED_KEYS if key in data})


def problem_args(problem: dict) -> list[str]:
    """Solver flags for a resolved "problem" object."""
    unknown = set(problem) - set(PROBLEM_FLAGS) - set(SWITCH_FLAGS) - {"blob_centre"}
    if unknown:
        raise ValueError(f"unknown problem keys: {sorted(unknown)}")
    args: list[str] = []
    for key, flag in PROBLEM_FLAGS.items():
        value = problem.get(key)
        if value is not None:
            args += [flag, format_value(value)]
    for key, flag in SWITCH_FLAGS.items():
        value = problem.get(key)
        if value is not None and not isinstance(value, bool):
            raise ValueError(f"{key} must be true or false")
        if value:
            args.append(flag)
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
