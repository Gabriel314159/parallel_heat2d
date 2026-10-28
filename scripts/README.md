# Experiment and plotting tools

| File | Status | Responsibility |
| --- | --- | --- |
| `make_decomposition_figures.py` | Implemented | Draw the row and 2D-block decomposition figures of `docs/design.md` into `docs/diagrams/` (Matplotlib only). |
| `check_correctness.py` | Implemented | Resolve the runs of `experiments/correctness.json`, run `heat_serial` on each and check mass conservation (relative mass change within `mass_conservation_rtol`), constant-field invariance under the flow, and `cosine_decay` against the exact discrete solution $(1-\lambda\Delta t)^n$; the error against the continuous solution is reported separately. With `heat_parallel` built, also run each case's parallel configuration and compare its field cell by cell with the serial reference (`--serial-only` skips this). Exits nonzero on a failed check. |
| `plot_field.py` | Implemented | Plot one field (a final pigment field, or one snapshot: pigment, vorticity or $\psi$ with `--field`) on its physical coordinates, with optional streamlines or arrows of the flow from the snapshot's $\psi$ (`--flow`). |
| `animate_field.py` | Implemented | Turn a snapshot folder into a GIF or MP4 of the pigment, the vorticity or $\psi$ (`--field`), with a fixed colour scale, the time and the mass change or energy in the title, and optional streamlines (`--flow`). |
| `heat2d_tools.py` | Implemented | Shared helpers: merge the JSON configuration, build solver command lines, parse the solver's `RESULT` line. |
| `fieldviz.py` | Implemented | Shared plotting helpers: read the output formats, colour scales, streamlines and arrows from $\psi$. |
| `run_benchmarks.py` | Planned | Expand explicit experiment resources, preserve resolved parameters, run warmups/repetitions and append raw measurement records. |
| `plot_results.py` | Planned | Read measured CSV files and produce timing, speedup, efficiency and strategy-comparison figures. |

## Requirements

- `check_correctness.py` uses only the Python standard library (Python 3.9+), so it runs wherever the solver does. CTest runs it as the `correctness_serial` and `correctness_parallel` tests.
- The plotting scripts need NumPy and Matplotlib 3.6+: `python3 -m pip install -r requirements.txt`. MP4 output also needs `ffmpeg`.
- Build the solver first: `cmake -S . -B build && cmake --build build`. The scripts look for it in `build/` (option `--build-dir`).

## Usage

```bash
# All correctness checks, serial and parallel (fields go to outputs/correctness/)
python3 scripts/check_correctness.py
python3 scripts/check_correctness.py --launcher "mpiexec --oversubscribe -n {ranks}"   # fewer than 6 cores

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

## Runner requirements (for the planned benchmark tools)

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
