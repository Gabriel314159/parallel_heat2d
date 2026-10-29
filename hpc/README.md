# Atlas / Vanda run setup

This project has not yet been configured or measured on either machine. The assignment requires runtime figures from Atlas or Vanda.

## Fill in after checking the current platform documentation

- Host or partition used, node CPU model and physical cores per node.
- Scheduler, queue, project/account options and requested wall time.
- Compiler, MPI and OpenMP module versions.
- Approved MPI launcher and rank/thread placement syntax.
- Number of nodes, ranks per node, threads per rank and memory requested.
- Appropriate filesystem for builds, temporary outputs and retained measurements.

Do not assume that a generic `mpiexec` command or scheduler header is valid on either platform. The supplied [job_template.sh](job_template.sh) has no scheduler-specific directives or module names. Add verified directives before submitting it.

## Large single runs

[big_run.sh](big_run.sh) runs one large case with snapshots, with `heat_serial` or with `heat_parallel` on `RANKS` × `THREADS` cores; its settings come from environment variables. Add the scheduler directives for the platform and set `LAUNCHER` (for example `srun`) before submitting it on a cluster. The main README, section "Larger runs on a bigger machine", explains how the cost grows with the grid. Such runs are for looking at the dynamics; the timed measurements for the report follow the run sequence below.

## Run sequence after implementation

1. Build in Release mode with the selected compiler and MPI environment.
2. Run the correctness checks (experiment 0, `scripts/check_correctness.py`) inside an allocation. The timings are only believed if they pass with this build on this machine.
3. Verify the CPU allocation, actual thread count and rank placement. The run report and the `RESULT` line give the ranks, threads and nodes the solver actually used; binding is not reported and must be checked separately.
4. Run experiments 1 to 5 of [experiments/README.md](../experiments/README.md), in that order, with the planned warmups and measured repetitions and without unrelated concurrent work. Time one run of the benchmark problem first, and adjust its size in `experiments/benchmark_problem.json` if needed.
5. Copy compact raw metrics and environment records into `results/raw/`.
6. Generate figures from those metrics and cite the exact commit/configuration in the report.

For the first hybrid implementation, test thread binding explicitly. Settings such as `OMP_PLACES=cores` and `OMP_PROC_BIND=close` are candidates to verify against the actual node architecture, not universal optimum settings.

Record the scheduler job ID, launcher, module list and hardware details with the results. Large solver snapshots belong in temporary output storage, not Git.
