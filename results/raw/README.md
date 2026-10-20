# Raw measurements

No performance measurements have been collected yet.

[metrics_columns.csv](metrics_columns.csv) is a header-only format template. It must not be treated as experimental data.

Commit compact raw CSV files and associated resolved configurations/environment records when experiments are run. Keep every measured repetition, its success/failure status and the exact Git commit. A row's `elapsed_s` is the maximum solver-loop time across ranks for a parallel run, or the serial solver-loop time for the serial reference.

`total_cores = mpi_ranks * omp_threads`; record `nodes`, placement and the actual observed thread count in session metadata. `max_abs_error` is optional for performance-only runs and should be blank when not checked. Analytical discretization errors need a separate validation record.

Do not commit full-resolution temperature dumps, build directories or temporary log bundles here. Large temporary files belong in ignored output directories or the appropriate HPC filesystem.
