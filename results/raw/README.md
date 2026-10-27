# Raw measurements

No performance measurements have been collected yet.

[metrics_columns.csv](metrics_columns.csv) is a header-only format template. It must not be treated as experimental data.

Commit compact raw CSV files and associated resolved configurations/environment records when experiments are run. Keep every measured repetition, its success/failure status and the exact Git commit. A row's `elapsed_s` is the maximum solver-loop time across ranks for a parallel run, or the serial solver-loop time for the serial reference.

`total_cores = mpi_ranks * omp_threads`; record `nodes`, placement and the actual observed thread count in session metadata. `max_abs_error` is optional for performance-only runs and should be blank when not checked. `mass_error` is the change of the total scalar relative to its initial absolute amount, $\lvert M_{\text{end}}-M_0\rvert / (\Delta x\,\Delta y\sum\lvert c^0_{i,j}\rvert)$, computed outside the timed loop. For non-negative fields this is $\lvert M_{\text{end}}-M_0\rvert/M_0$; unlike that ratio, it is also defined for `cosine_decay`, where $M_0 = 0$. `peclet` and `peclet_grid` are $\mathrm{Pe}$ and $\mathrm{Pe}_{\text{grid}}$ as defined in `docs/design.md`. Analytical discretization errors need a separate validation record.

Do not commit full-resolution field dumps, build directories or temporary log bundles here. Large temporary files belong in ignored output directories or the appropriate HPC filesystem.
