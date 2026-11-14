# Experiment and plotting tools

| File | Status | Responsibility |
| --- | --- | --- |
| `make_decomposition_figures.py` | Implemented | Draw the row and 2D-block decomposition figures of `docs/design.md` into `docs/diagrams/` (Matplotlib only). |
| `check_correctness.py` | Implemented | Resolve the runs of `experiments/correctness.json`, run `heat_serial` on each and check mass conservation (relative mass change within `mass_conservation_rtol`), constant-field invariance under the flow, and `cosine_decay` against the exact discrete solution $(1-\lambda\Delta t)^n$; the error against the continuous solution is reported separately. With `heat_parallel` built, also run each case's parallel configuration and compare its field cell by cell with the serial reference (`--serial-only` skips this). Exits nonzero on a failed check. |
| `plot_field.py` | Implemented | Plot one field (a final pigment field, or one snapshot: pigment, vorticity or $\psi$ with `--field`) on its physical coordinates, with optional streamlines or arrows of the flow from the snapshot's $\psi$ (`--flow`). |
| `animate_field.py` | Implemented | Turn a snapshot folder into a GIF or MP4 of the pigment, the vorticity or $\psi$ (`--field`), with a fixed colour scale, the time and the mass change or energy in the title, and optional streamlines (`--flow`). |
| `heat2d_tools.py` | Implemented | Shared helpers: merge the JSON configuration (an experiment inherits its `defaults` file, which may inherit another), build solver command lines, parse the solver's `RESULT` line. |
| `fieldviz.py` | Implemented | Shared plotting helpers: read the output formats, colour scales, streamlines and arrows from $\psi$. |
| `run_benchmarks.py` | Implemented | Run one performance experiment (`--config experiments/<name>.json`), adapted from the main branch's runner. It expands the file as `experiments/README.md` describes ("Run expansion"), runs every configuration `warmup_runs` + `repetitions` times one after another, and appends one row per run, warmups included, to the file's `metrics_csv`. It writes the resolved plan, the environment and every run's output to `<output.directory>/<session>/`. Failed runs are kept as rows with status `error`. |
| `plot_results.py` | Planned | Read measured CSV files and produce timing, speedup, efficiency and strategy-comparison figures. |

## Requirements

- `check_correctness.py` and `run_benchmarks.py` use only the Python standard library (Python 3.9+), so they run wherever the solver does. CTest runs `check_correctness.py` as the `correctness_serial` and `correctness_parallel` tests, and on `experiments/correctness_3d.json` as `correctness_3d_serial` and `correctness_3d_parallel`. An experiment file with `"model": "3d"` selects the 3D programs in both scripts (`MODELS` in `heat2d_tools.py`). `plot_field.py` also draws the vertical sections of 3D snapshots (`section_*.npy`).
- The plotting scripts need NumPy and Matplotlib 3.6+: `python3 -m pip install -r requirements.txt`. MP4 output also needs `ffmpeg`.
- Build the solver first: `cmake -S . -B build && cmake --build build`. The scripts look for it in `build/` (option `--build-dir`).

## Usage

```bash
# All correctness checks, serial and parallel (fields go to outputs/correctness/)
python3 scripts/check_correctness.py
python3 scripts/check_correctness.py --launcher "mpiexec --oversubscribe -n {ranks}"   # fewer than 6 cores

# A performance experiment: first the plan (configurations and commands), then the runs
python3 scripts/run_benchmarks.py --config experiments/strong_scaling_mpi.json --dry-run

# The 3D pigment model: the files with "model": "3d" run heat3d_serial and heat3d_parallel
python3 scripts/check_correctness.py --experiment experiments/correctness_3d.json
python3 scripts/run_benchmarks.py --config experiments/strong_scaling_3d.json --dry-run
python3 scripts/run_benchmarks.py --config experiments/strong_scaling_mpi.json --run serial_reference --run p4_t1_2d_blocking

# Snapshots: pigment, psi and vorticity every 500 steps
OMP_NUM_THREADS=2 mpiexec -n 4 ./build/heat_parallel --decomp 2d \
    --snapshot-interval 500 --snapshot-coarsen 2 --snapshot-dir outputs/preview

# One snapshot: the pigment with streamlines, or the vorticity
python3 scripts/plot_field.py outputs/preview/field_025000.npy --flow streamlines --gamma 0.5
python3 scripts/plot_field.py outputs/preview/field_025000.npy --field vorticity

# Animations
python3 scripts/animate_field.py outputs/preview --flow streamlines -o outputs/preview.gif
python3 scripts/animate_field.py outputs/preview --field vorticity -o outputs/vorticity.mp4
```

Every script has `--help`. Useful options:

- `plot_field.py --gamma 0.5` and `animate_field.py --gamma G`: a power colour scale that keeps faint, spread-out pigment visible (the animation uses 0.5 by default; `--gamma 1` is linear). The colour bar states the scale.
- `animate_field.py --every K` uses every K-th snapshot; `--frames-dir DIR` also saves each frame as PNG.
- `check_correctness.py --run NAME` runs selected cases only; `--launcher` sets the MPI launch command for the parallel comparisons, for example `--launcher "srun -n {ranks}"`.

Colours: non-negative fields use one blue scale from near-white (no pigment) to dark blue; signed fields such as the vorticity, $\psi$ or `cosine_decay` use a blue/red scale with a neutral grey at zero. The animation keeps one colour scale for all frames: for the pigment, fading by diffusion stays visible, and by the discrete maximum principle no frame leaves the range of the first; for the vorticity and $\psi$, the scale covers the largest value of all frames. Streamlines are drawn as contours of $\psi$, which are exactly the streamlines of the flow.

`run_benchmarks.py` options for a cluster:

- `--launcher` and `--mpi-arg=VALUE` (repeated) build the launch command, as on the main branch. The options may contain placeholders: `{ranks}` and `{threads}`; with `--cores-per-node`, `{nodes}` and `{ranks_per_node}`; with `--sockets-per-node` as well, `{ranks_per_socket}`. A placeholder that is not a whole number for some configuration stops the plan before anything runs.
- `--launch-serial` starts the serial reference through the launcher with one task, so that it is bound to a core like the other runs.
- `--require-binding` counts a run as failed unless it reports `binding` `cores`.
- `--allocated-cores N` declares the cores of a multi-node allocation (by default, the cores this process may use).
- `--run NAME` (repeated) runs only some configurations, for example to split an experiment over jobs; `--dry-run` lists the names. `--keep-going` continues after a failed run; `--timeout` limits each run (default 3600 s).

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

The solver's `RESULT` line already contains most metrics columns. The expected raw CSV columns are listed in [results/raw/metrics_columns.csv](../results/raw/metrics_columns.csv). The file contains a header only and is not an experimental result.
