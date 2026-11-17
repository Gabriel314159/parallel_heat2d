# Model B (3D flow with pressure projection), 1024^3 strong scaling

Source: the Model B package (`model_b.cpp`, 3D MAC grid, AB3 predictor, Neumann pressure Poisson by CG at every step), not in this repository. Run on Vanda on 16-17 Nov 2026, one thread per rank, unit cube of 1024^3 cells, one time step from rest, CG tolerance 1e-8, `--map-by core` up to 32 ranks and 32 ranks per socket at 64. The files in `jobs/` are the scripts used. `validation_*.csv` are the 69 correctness checks each job passed before timing anything.

| Ranks | Job | Time, s | Speedup | Efficiency | ms per CG it |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | 1567816 | 41151 (extrapolated) | 1 | 1 | 11641 |
| 2 | 1567815 | 20760 | 1.98 | 0.99 | 5873 |
| 4 | 1567815 | 10470 | 3.93 | 0.98 | 2962 |
| 8 | 1567814 | 5476 | 7.51 | 0.94 | 1549 |
| 16 | 1567814 | 3235 | 12.72 | 0.80 | 915 |
| 32 | 1567814 | 2350 | 17.51 | 0.55 | 665 |
| 64 | 1567814 | 1173 (median of 3: 1170-1176) | 35.08 | 0.55 | 332 |

Every run takes 3535 CG iterations, reaches a relative residual of 9.88e-9, a divergence of 1.2e-11 and a mass error below 5e-14.

**The one-core run did not finish.** The runner's 41400 s limit stopped it at iteration 3500 of 3535, 35 iterations short, so there is no `raw.csv` for `strong_p1`, only its log. The time above is the log's elapsed time at iteration 3500 (40746 s) plus 35 iterations at the rate of the last 200 (11.59 s), so it is good to about 1%, but it is not a measurement. Speedup and efficiency use it.

**Placement.** Up to 32 ranks all sit on socket 0 (`--map-by core`), as in the 2D experiments, so one socket's memory bandwidth is shared: the efficiency falls from 0.94 at 8 ranks to 0.55 at 32, then the step from 32 to 64 ranks, which uses both sockets, is 2.0x. The two-socket cost of 16 and 32 ranks was not measured here.

**Not run yet:** 128 and 256 ranks and the weak-scaling ladder (job 1567813, four nodes, still queued).
