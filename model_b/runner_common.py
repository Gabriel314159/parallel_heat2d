"""Common subprocess and result handling. Standard library only, Python >=3.8."""
import csv
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import time


def options(values):
    return [part for key, value in values.items() for part in ("--" + key.replace("_", "-"), str(value))]


def run_case(build, launcher, output, name, values, ranks=0, threads=1, timeout=600):
    build, output = Path(build).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    command = [str(build / ("model_b_mpi" if ranks else "model_b_serial"))] + options(values)
    if ranks:
        command = shlex.split(launcher) + ["-np", str(ranks)] + command
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE")
    # Serial binary has no OpenMP dependency and is genuinely one thread.
    started = time.perf_counter()
    with (output / (name + ".log")).open("w") as stream:
        stream.write("COMMAND " + shlex.join(command) + "\n")
        stream.flush()
        proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        try:
            proc.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
            # Stop only the subprocess group created for this run, including its
            # local MPI children. No user's other processes/jobs are touched.
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            raise RuntimeError("interrupted or timed out: " + name + "; see " + str(output / (name + ".log"))) from exc
    end_to_end = time.perf_counter() - started
    if proc.returncode:
        raise RuntimeError("failed: " + name + "; see " + str(output / (name + ".log")))
    lines = (output / (name + ".log")).read_text().splitlines()
    results = [json.loads(line[7:]) for line in lines if line.startswith("RESULT ")]
    if len(results) != 1:
        raise RuntimeError("expected exactly one RESULT in " + name)
    result = results[0]
    if result["ranks"] != max(1, ranks) or result["threads"] != threads:
        raise RuntimeError("actual MPI/OpenMP layout differs from request: " + name)
    result.update(case=name, end_to_end_s=end_to_end)
    return result


def read_field(path):
    with Path(path).open() as f:
        rows = list(csv.DictReader(f))
    data = {}
    for row in rows:
        key = tuple(float(row[k]) for k in ("x", "y", "z"))
        if key in data:
            raise ValueError("duplicate cell coordinate in " + str(path))
        data[key] = {k: float(row[k]) for k in ("c", "phi", "u_e", "v_n", "w_t")}
    return data


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
