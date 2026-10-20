# Experiment tools to implement

The scripts listed below are planned. No Python implementations are included yet.

| Planned file | Responsibility |
| --- | --- |
| `check_correctness.py` | Resolve correctness cases, invoke matching serial/parallel runs, compare full fields and exit nonzero on failed checks. |
| `run_benchmarks.py` | Expand explicit experiment resources, preserve resolved parameters, run warmups/repetitions and append raw measurement records. |
| `plot_results.py` | Read measured CSV files and produce timing, speedup, efficiency and strategy-comparison figures. |
| `plot_temperature.py` | Plot a field matrix and its physical coordinates, using comparable temperature scales. |

## Runner requirements

- Use `subprocess.run` with an argument list, not shell-interpolated command text.
- Set thread environment variables explicitly and preserve launcher/resource metadata.
- Run locally only within available CPU and memory resources. On HPC, run within the scheduler allocation.
- Record failed runs as failures; do not synthesize zero timings or silently drop them.
- Save all raw repetitions, plus a resolved configuration and environment metadata per experiment session.
- Keep matrix output and verification outside measured solver-loop time.
- Do not run all strong-scaling configurations concurrently on the same allocated resources.
- For overlap experiments, do not add independently measured communication and computation times to claim total time: those phases can overlap.
- Keep analytical-discretization error separate from serial/parallel equivalence error.

The expected raw CSV columns are listed in [results/raw/metrics_columns.csv](../results/raw/metrics_columns.csv). The file contains a header only and is not an experimental result.
