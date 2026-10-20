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

## Run sequence after implementation

1. Build in Release mode with the selected compiler and MPI environment.
2. Run a small serial/parallel validation case inside an allocation.
3. Verify the CPU allocation, actual thread count and rank placement.
4. Run planned warmups and measured repetitions without unrelated concurrent work.
5. Copy compact raw metrics and environment records into `results/raw/`.
6. Generate figures from those metrics and cite the exact commit/configuration in the report.

For the first hybrid implementation, test thread binding explicitly. Settings such as `OMP_PLACES=cores` and `OMP_PROC_BIND=close` are candidates to verify against the actual node architecture, not universal optimum settings.

Record the scheduler job ID, launcher, module list and hardware details with the results. Large solver snapshots belong in temporary output storage, not Git.
