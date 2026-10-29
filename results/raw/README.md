# Raw measurements

No performance measurements have been collected yet.

[metrics_columns.csv](metrics_columns.csv) is a header-only format template. It must not be treated as experimental data.

Commit compact raw CSV files and associated resolved configurations/environment records when experiments are run. Each experiment of `experiments/` has its own file, named after its experiment file (for example `strong_scaling_mpi.csv`). Keep every measured repetition, its success/failure status and the exact Git commit. A row's `elapsed_s` is the maximum solver-loop time across ranks for a parallel run, or the serial solver-loop time for the serial reference.

`total_cores = mpi_ranks * omp_threads`; `nodes` and the actual observed thread count (`omp_threads`) come from the solver's `RESULT` line; record the placement in session metadata. `cores_per_rank` is the fewest physical cores any rank's threads may run on (0 if the system does not report it); a row with `cores_per_rank` below `omp_threads` had threads sharing cores and is not a valid timing. `poisson_s` is the time of the Poisson solves in the timed loop; `poisson_halo_s` and `poisson_sum_s` are the parts of it spent inside the CG's halo exchanges (under overlap, only posting and waiting) and global sums. All three are maxima over the ranks, so they do not add up. `max_abs_error` is optional for performance-only runs and should be blank when not checked. `mass_error` is the change of the total scalar relative to its initial absolute amount, $\lvert M_{\text{end}}-M_0\rvert / (\Delta x\,\Delta y\sum\lvert c^0_{i,j}\rvert)$, computed outside the timed loop. For non-negative fields this is $\lvert M_{\text{end}}-M_0\rvert/M_0$; unlike that ratio, it is also defined for `cosine_decay`, where $M_0 = 0$. `peclet` and `peclet_grid` are $\mathrm{Pe}$ and $\mathrm{Pe}_{\text{grid}}$ as defined in `docs/design.md`. Analytical discretization errors need a separate validation record.

Do not commit full-resolution field dumps, build directories or temporary log bundles here. Large temporary files belong in ignored output directories or the appropriate HPC filesystem.
