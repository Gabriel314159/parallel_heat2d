# Four-node MPI audit and diagnostic update

This directory contains the MPI runner and node probe used by `four_nodes_checked.pbs`.

## Run on Vanda

Submit from `poisson3d_benchmark/` with the name of your Vanda allocation:

```bash
qsub -P YOUR_PROJECT four_nodes_checked.pbs
```

The job uses the original build.sh and poisson3d.cpp. A fresh build and output
directory are used for each PBS job. One warm-up and five measured repetitions,
all problem sizes, solver tolerances and CSV columns remain unchanged.

## Findings and scope

- Allocation is four distinct, exclusive nodes, 72 CPU slots each. The 256-core
  comparisons use 64 active cores per node. This is not oversubscription.
- The original `--prefix EBROOTPRRTE`, library forwarding, MPI/OpenMP environment
  forwarding and six socket-based hybrid mappings match the 2D implementation.
- The 4-rank/64-thread case uses one rank per node and explicit OpenMP places
  spanning both sockets, as in 2D. It requires the documented CPU numbering.
- Correct weak scaling at 36 ranks to 18 ranks per socket and at 144 ranks to
  24 ranks per socket on three nodes. The original fixed 32-per-socket maximum
  did not enforce these balanced layouts. Other mappings are unchanged.
- This correction is not a demonstrated fix for the reported PRRTE lost-daemon
  failure. The supplied log does not identify its root cause.
- Do not change MPI implementation or network-interface flags without evidence.
  The PBS-node file is read for inventory only and is not modified.

Reference code:
https://github.com/Gabriel314159/parallel_heat2d/blob/43ad716291c270eb6cfbbb2d99ccd8a374f8fcc4/hpc/jobs/common.sh
https://github.com/Gabriel314159/parallel_heat2d/blob/43ad716291c270eb6cfbbb2d99ccd8a374f8fcc4/hpc/jobs/mpirun.sh
https://github.com/Gabriel314159/parallel_heat2d/blob/6f53d23b06227deabef691b31bd4edbc7308edcb/experiments/weak_scaling_extra.json

## Same-job checks and output

Before correctness validation and measurements, the runner launches one probe
per allocated node. It checks the CPU numbering used by the 64-thread mapping.
Then it solves a small manufactured Poisson problem with one MPI rank per node,
checking four distinct nodes, thread binding and numerical error. Each of these
two launches has a maximum timeout of 120 seconds. This does not predict queue
waiting time or guarantee that later high-rank jobs will succeed.

If the checks pass, all four-node experiments run immediately in the same PBS
allocation. If a check fails, the job stops with diagnostic logs and does not
publish that check as a timing sample. There is no automatic retry or checkpoint
resume in this update.

Files under `results/<PBS_JOBID>/`:

- `allocation.json`: allocated hosts, slots and modules, written before validation.
- `launch-logs/launch-NNNN.json`: exact command and relevant environment.
- `launch-logs/launch-NNNN.stdout` and `.stderr`: completed launch output, including
  MPI binding reports. Captured partial output is saved on a launcher timeout.
- `preflight-nodes.json`, `preflight-solver.json`: successful preflight evidence.
- `four-nodes.json`, `four-nodes.csv`: original metadata/results format, created
  after the correctness gate. Each completed sample is flushed to the CSV.

Do not combine incomplete or failed measurements as though the campaign finished.
The existing summarize.py remains compatible. Original one-node and GPU results
do not need to be rerun because of this launch-script update.

## Local verification

Python compilation and PBS shell syntax checks passed. All 76 generated layouts
(including optional controls and CPU/GPU suites) were compared against the original
runner: only weak36/144 mappings change, plus binding-report diagnostics. Resource
counts and specified process grids were checked. Log preservation was exercised
for successful launches, nonzero exits, and simulated timeouts. Four-node runtime
and Linux topology checks have not been executed on Vanda by this assistant.
