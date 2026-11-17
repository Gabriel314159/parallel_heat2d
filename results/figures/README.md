# Report figures

Made by `scripts/plot_*.py` from the CSVs of `results/raw/`; every figure is saved as `.png` (200 dpi) and `.pdf` here. Run a script from `scripts/`, for example `python3 plot_strong_mpi.py`.

- `strong_scaling_mpi.png`: experiment 1, runtime, speedup and efficiency against `heat_serial`.
- `strong_scaling_mpi_iteration.png`: experiment 1, time per CG iteration.
- `strong_scaling_openmp.png`: experiment 2, runtime, speedup and efficiency against `heat_serial`.
- `strong_scaling_openmp_iteration.png`: experiment 2, time per CG iteration.
- `weak_scaling_mpi.png`: experiment 5, runtime, time per CG iteration and the two weak efficiencies, from `weak_scaling_rerun.csv`.
- `weak_scaling_mpi_scaled.png`: experiment 5 against the serial program on the same grid: runtime, scaled speedup $S_s=T_{\text{serial}}(N_p)/T(p)$ and $E_s=S_s/p$ with $E_{w,\text{it}}$, from `weak_scaling_rerun.csv` and `serial_weak_*.csv`, by `plot_weak_scaled.py`.
- `strong_scaling_3d.png` and `strong_scaling_3d_iteration.png`: Model B (3D), 1024³ cells, runtime, speedup, efficiency and time per CG iteration, from `model_b_3d/strong_p*/`, by `plot_strong_3d.py`; the 1-rank point is extrapolated.
- `domain.png`: the basin, its walls, the initial pigment blob and the wind forcing against $y$ (report Section 2), by `plot_domain.py` from `config/defaults.json`.
- `domain_3d.png`: the 3D pigment as clouds of points at τ = 0, 2.44 and 9.77 (depth stretched 6 times), by `plot_domain_3d.py` from the fields in `outputs/flow3d` (heat3d_parallel of the 3D-Poisson branch, 128 × 128 × 8 cells, lz = 0.0625, viscosity 0.004, spin-up 8; the commands are in the .log files there).
- `flow_evolution.png`: the flow at release and the pigment at three times after it (report Section 4.2), by `plot_flow_evolution.py` from the snapshots in `outputs/flow128` (128 × 128 cells, viscosity 0.004, spin-up 8, 16000 steps, a snapshot every 1000 steps).
- `hybrid_layouts_iteration.png`: experiment 3, the time per CG iteration of the seven layouts and its two parts, the global sums and the rest. Its y axis is linear, not log (rule 4), so that the two parts visibly add up to $t_{\text{it}}$.

## Rules for every plot

The shared code is `scripts/plot_style.py`; a new plot imports it and follows these rules.

1. **Curves only.** Lines with markers, one panel per quantity. No bars, no dual axes.
2. **One entity, one colour, in every figure.** Blue is the measured curve (the baseline in speedup plots), orange the second serial reference or the second variant, aqua halo, yellow global sums. A new variant takes the next slot of `ROLE` (never recycled). The ideal curve (perfect scaling, labelled with its formula) is always dashed grey with no marker; a second reference (for example the runtime if only the CG iterations grow) is dotted grey.
3. **x axis.** Cores or ranks on a log2 axis with the plain integers as ticks (`ranks_axis`).
4. **y axis.** Log scale for times and speedups, linear 0 to 1.25 for efficiency, with the ideal line at 1. The label carries the unit: `runtime $T$ (s)`.
5. **Quantities are the report's:** $T$ is the median `elapsed_s` of the 5 measured repetitions, with a band from the minimum to the maximum; $S(p)=T_{\text{serial}}/T(p)$, $E(p)=S(p)/p$ against `heat_serial` (the baseline is named in the title); $t_{\text{it}}=T_P/(n\bar k)$ with the halo and sum parts. Warmup rows are never plotted.
6. **Text.** The title names the experiment as in the report ("Strong scaling with MPI"); in a multi-panel figure each panel has a one-word title (Runtime, Speedup, Efficiency). No note below the figure: the problem, the steps, the baseline and the statistic go in the report's caption. A legend whenever there are two or more curves.
7. **Size.** 11.5 × 3.6 in for three panels in a row, 5.5 × 3.8 in for one panel. Use `layout="constrained"`.
8. **Data.** Read only `results/raw/*.csv` (`ps.load`); join files by column name; never edit a CSV for a figure.
