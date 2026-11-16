#!/usr/bin/env python3
"""Summarize measured rows from the Vanda CPU and GPU suites."""

import argparse
import csv
import statistics
import sys
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", help="one-node.csv, four-nodes.csv, and/or gpu.csv from one experiment campaign")
    args = parser.parse_args()
    keys = ("experiment", "backend", "nx", "ny", "nz", "lx", "ly", "lz",
            "ranks", "threads", "nodes", "px", "py", "exchange", "gpu_count", "N", "nnz", "tol")
    metrics = ("solve_s", "solve_total_s", "device_setup_s", "device_finish_s",
               "halo_s", "halo_device_copy_s", "reduction_s", "reduction_min_s",
               "iterations", "solve_s_per_iteration")
    groups = defaultdict(list)
    for filename in args.csv:
        with open(filename, newline="") as handle:
            for row in csv.DictReader(handle):
                if row["phase"] == "measured" and row["binding"] == "cores":
                    if row["backend"] == "gpu" and row.get("poisson_device") != "gpu":
                        parser.error(f"{filename} includes a GPU row without confirmed device execution")
                    # CPU files from before offloading only recorded solve_s.
                    row.setdefault("solve_total_s", row["solve_s"])
                    row.setdefault("device_setup_s", "0")
                    row.setdefault("device_finish_s", "0")
                    row.setdefault("exchange", "blocking")
                    row.setdefault("gpu_count", "1" if row["backend"] == "gpu" else "0")
                    row.setdefault("halo_device_copy_s", "0")
                    row.setdefault("reduction_min_s", "")
                    groups[tuple(row[k] for k in keys)].append(row)
    if not groups:
        parser.error("no measured rows with valid binding")
    summary = []
    for key, runs in groups.items():
        row = dict(zip(keys, key))
        row["samples"] = len(runs)
        row["cores"] = int(row["ranks"]) * int(row["threads"])
        for metric in metrics:
            values = [float(run[metric]) for run in runs if run.get(metric)]
            row[metric] = statistics.median(values) if values else ""
        row["solve_min_s"] = min(float(run["solve_s"]) for run in runs)
        row["solve_max_s"] = max(float(run["solve_s"]) for run in runs)
        summary.append(row)

    strong = next((r["solve_s"] for r in summary if r["experiment"] == "strong" and r["backend"] == "serial"), None)
    openmp = next((r["solve_s"] for r in summary if r["experiment"] == "openmp" and r["backend"] == "serial"), None)
    weak = next((r["solve_s"] for r in summary if r["experiment"] == "weak" and r["backend"] == "mpi" and r["ranks"] == "1"), None)
    hybrid = {r["experiment"]: r["solve_s"] for r in summary
              if r["experiment"] in ("hybrid", "hybrid-node") and r["threads"] == "1"}
    problem_keys = ("nx", "ny", "nz", "lx", "ly", "lz", "tol")
    cpu_gpu_references = {
        tuple(r[k] for k in problem_keys + ("ranks", "threads", "px", "py", "exchange")): r["solve_total_s"]
        for r in summary if r["experiment"] == "gpu-size" and r["backend"] == "cpu-omp"
    }
    one_gpu = {tuple(r[k] for k in problem_keys + ("threads",)): r["solve_total_s"]
               for r in summary if r["experiment"] == "gpu-size" and r["gpu_count"] == "1"}
    for r in summary:
        r["speedup"] = r["efficiency"] = r["speedup_vs_pure_mpi"] = ""
        r["gpu_speedup_vs_matched_cpu"] = r["gpu_speedup_vs_one_gpu"] = ""
        if r["experiment"] == "strong" and strong is not None:
            r["speedup"] = strong / r["solve_s"]
            r["efficiency"] = r["speedup"] / r["cores"]
        elif r["experiment"] == "openmp" and openmp is not None:
            r["speedup"] = openmp / r["solve_s"]
            r["efficiency"] = r["speedup"] / r["cores"]
        elif r["experiment"] == "weak" and weak is not None:
            r["efficiency"] = weak / r["solve_s"]
            r["speedup"] = r["cores"] * r["efficiency"]
        elif r["experiment"] in hybrid:
            r["speedup_vs_pure_mpi"] = hybrid[r["experiment"]] / r["solve_s"]
        elif r["experiment"] == "gpu-size" and r["backend"] == "gpu":
            match = tuple(r[k] for k in problem_keys + ("ranks", "threads", "px", "py", "exchange"))
            if match in cpu_gpu_references:
                r["gpu_speedup_vs_matched_cpu"] = cpu_gpu_references[match] / r["solve_total_s"]
            reference = tuple(r[k] for k in problem_keys + ("threads",))
            if reference in one_gpu:
                r["gpu_speedup_vs_one_gpu"] = one_gpu[reference] / r["solve_total_s"]
    writer = csv.DictWriter(sys.stdout, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)


if __name__ == "__main__":
    main()
