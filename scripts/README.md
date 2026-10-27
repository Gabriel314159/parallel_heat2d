# Experiment tools

All listed Python tools are implemented. Run commands from the repository root with the Python environment activated.

| File | Responsibility |
| --- | --- |
| `common.py` | Recursive configuration merging, stable-step checks, MPI launching, JSON parsing and provenance. |
| `check_correctness.py` | Hand-computed update, physical invariants, serial/parallel field comparison and sine refinement. |
| `run_benchmarks.py` | Explicit experiment expansion, resource guard, warmups, raw repetitions and logs. |
| `plot_temperature.py` | A physical-coordinate temperature PNG from a final CSV field. |
| `plot_results.py` | Latest-session measured medians, min–max runtime ranges, speedup and efficiency. |

Use `--help` on each command. MPI tools accept `--build-dir`, `--launcher`, repeated `--mpi-arg=VALUE` and `--timeout`. Launcher arguments expand `{threads}` and `{ranks}` per run. This accommodates platform-specific placement without shell interpolation.

The benchmark runner checks requested workers against local CPU affinity unless `--allocated-cores` explicitly declares an allocation. It does not allocate cluster resources. Warmups and failed runs are retained with distinct phase/status fields; failed solver invocations exit nonzero. The plotter excludes them from measured statistics and keeps sessions separate.

Resolved configurations and environment records are written beside per-run stdout/stderr files. The CSV format is listed in [metrics_columns.csv](../results/raw/metrics_columns.csv). That file is a header template, not benchmark data.

Only final fields are supported; `snapshot_interval` must remain zero. Configure fixed boundaries explicitly when changing JSON cases from the heated-top defaults. Analytical error is separate from serial/parallel equivalence error. Phase times are not added to claim total time under overlap.
