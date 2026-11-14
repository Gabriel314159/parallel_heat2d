# Atlas / Vanda run setup

This project has not yet been measured on either machine. The assignment requires runtime figures from Atlas or Vanda; we use Vanda.

## Vanda

[check_platform.sh](check_platform.sh) reports a machine's nodes, queues and modules and changes nothing. It ran on Vanda's login node on 3 November 2026, from the laptop:

```bash
ssh <user>@vanda.nus.edu.sg 'bash -l -s' < hpc/check_platform.sh > outputs/vanda_platform.txt
```

The login node is a virtual machine of another type than the compute nodes, and the site asks us to compile on a compute node. Inside a job, the same script reports the compute node and tests which compilers build OpenMP offload code there. It ran on 3 November 2026 in two short jobs of the free queue, one on a CPU node (`cn-161`, 2 cores) and one on a GPU node (`gn-a40-004`, 4 cores and 1 GPU); the reports are `outputs/vanda_platform_cpu.txt` and `outputs/vanda_platform_gpu.txt`. A job script must not give `#PBS -o` a relative path: PBS reads it from the submission directory, and a job whose output directory does not exist is restarted again and again (exit status −3) without running.

- **System:** `vanda.nus.edu.sg`, Red Hat Enterprise Linux 9.4, PBS Pro (`qsub`, `qstat`). `hpc cpu` and `hpc gpu` print sample job scripts.
- **CPU nodes** (`cn-…`, about 140): 72 cores and 501 GB each. `lscpu` on a CPU node: two Intel Xeon Platinum 8452Y of 36 cores, one thread per core (no hyper-threading), two NUMA nodes, one per socket. The core numbers alternate between the sockets: the even ones are socket 0, the odd ones socket 1 (see "Binding").
- **Memory pages:** transparent huge pages (2 MiB) are `always` on the login node; the log of the follow-up jobs ([jobs/extra_threads.pbs](jobs/extra_threads.pbs) and the others) prints the compute node's setting. `heat_parallel` turns them off for itself when it runs more than one thread, so that the first touch of a page by a thread places it on that thread's socket (`docs/design.md`, "Data storage and shared computation").
- **GPU nodes** (`gn-a40-…`, about 90): the same two CPUs, 501 GB and two NVIDIA A40 (46 GB usable, compute capability 8.6, driver 575.57). The site gives each GPU one socket of 36 cores; `nvidia-smi topo -m` puts GPU 0 on NUMA node 0.
- **OpenMP offload** (tested on the A40): `nvc++ -mp=gpu` from `nvhpc/26.1` and `clang++ -fopenmp --offload-arch=sm_80` from `LLVM/20.1.7-GCCcore-14.2.0` both build a test program whose target region runs on the GPU. Neither `g++` can: the system `g++` 11.4 names the target but lacks its offload tools, and `gcc/13.1.0` was built without it. On a CPU node the `nvc++` program runs its target region on the host, and `clang++` does not build (no CUDA libraries there).
- **Free queues** (`-q auto_free`): `cpu_parallel`, 2 to 160 cores, at most 2 running jobs; `gpu`, at most 36 cores and 2 GPUs, 48 hours.
- **Charged queues** (`-P <project>`): `batch_cpu`, up to 720 cores on 10 nodes, 240 hours, 6 running jobs; `batch_gpu`, up to 2 GPUs. The 256-core runs need `batch_cpu`, so a project with credits (`amgr login`, then `hpc project`).
- **Modules** (EasyBuild names): `GCC/14.3.0` with `OpenMPI/5.0.8-GCC-14.3.0`, `CMake/4.0.3-GCCcore-14.3.0`, `nvhpc/26.1` and `nvhpc-hpcx-cuda12/26.1` (NVIDIA's compilers, with HPC-X MPI), `CUDA/12.*`, `LLVM/20.1.7`, Intel compilers and Intel MPI.
- **Storage:** home 40 GB; `/scratch/$USER` 500 GB, where files idle for 60 days are deleted.

Do not assume that a generic `mpiexec` command or scheduler header is valid on either platform. The supplied [job_template.sh](job_template.sh) has no scheduler-specific directives or module names. Add verified directives before submitting it.

## Binding

Every timed run puts one thread on each physical core and leaves the other hardware threads of that core unused (no SMT, or hyper-threading). Two of our threads on one core share its arithmetic units and caches; that core would then count twice in the number of cores $p$, and the speedup and efficiency would be wrong.

1. **Count physical cores, not hardware threads.** `lscpu` shows `Thread(s) per core`, `Core(s) per socket` and `Socket(s)`; a node has sockets × cores per socket physical cores. Ranks per node × threads per rank must not exceed that. Where the scheduler counts hardware threads as CPUs, ask for whole cores.
2. **Bind each thread to a core.** The job scripts set
   ```bash
   export OMP_PLACES=cores
   export OMP_PROC_BIND=close
   ```
   Each OpenMP place is then one physical core, and the threads of a rank take one place each, on neighbouring cores. With `OMP_PLACES=threads`, close binding would put two threads on the two hardware threads of one core.

   On Vanda, neighbouring core numbers lie on different sockets. An OpenMP runtime that lists the places in the order of the core numbers would then spread 36 threads over both sockets, not fill one. Bind the rank first (step 3): its places are then only the cores the launcher gave it, and Open MPI counts cores socket by socket. Where one rank holds a whole node, list the places explicitly ("Placement on Vanda").
3. **Bind each rank to its cores.** Give every rank as many whole cores as it has threads, and bind it to them. The options differ between launchers; these are candidates, to check against the Atlas/Vanda documentation:
   - Slurm: `--ntasks=R --cpus-per-task=T --hint=nomultithread` for the job, and `srun --cpus-per-task=T --cpu-bind=cores` (recent Slurm versions do not pass `--cpus-per-task` on from the job to `srun`).
   - Open MPI's `mpirun`: `--map-by slot:PE=T --bind-to core`. Never `--use-hwthread-cpus`, which makes every hardware thread a slot.
   - PBS: `select=N:ncpus=C:mpiprocs=R:ompthreads=T`, with $C = R \times T$ physical cores per node, and the binding options of the MPI launcher as above.
4. **Bind the serial reference too**, to one core: start it through the launcher with one task (for example `srun -n 1 --cpu-bind=cores` or `mpirun -n 1 --bind-to core`).
5. **Check every run.** The solvers report `binding` in the run report and on the `RESULT` line:
   - `cores`: every thread is bound to a physical core of its own;
   - `shared`: two threads on one node are bound to the same core, for example to its two hardware threads, or a node runs more threads than it has cores;
   - `unpinned`: some thread may run on more than one core, so the operating system decides where it runs (usually `OMP_PROC_BIND` is unset, or a rank is not bound);
   - `unknown`: the system does not report the binding (macOS).

   Only runs with `binding` `cores` are valid timings; `heat_parallel` warns otherwise. For the first job of every layout, also look at the launcher's own report (`srun --cpu-bind=verbose,cores`, `mpirun --report-bindings`, or `OMP_DISPLAY_AFFINITY=TRUE`), and record the options used.

## Placement on Vanda

Decided on 3 November 2026; check it in the first job with `mpirun --report-bindings`, `OMP_DISPLAY_AFFINITY=TRUE` and the solver's `binding`, which must be `cores`.

- **Whole nodes.** Every timed run reserves its nodes exclusively (`#PBS -l place=excl`), so that no other job shares them. The project is charged for all 72 cores of each node.
- **Build.** `GCC/14.3.0`, `OpenMPI/5.0.8-GCC-14.3.0` and `CMake/4.0.3-GCCcore-14.3.0`, a Release build with `-DHEAT2D_NATIVE=ON`, compiled on a compute node. `heat_serial`, `heat_serial_omp` and `heat_parallel` come from the same build.
- **Core numbers alternate between the sockets.** The even core numbers 0, 2, …, 70 are socket 0 and the odd ones 1, 3, …, 71 socket 1. Open MPI's mappings count cores socket by socket (hwloc's logical order), so `--map-by core` and `ppr:N:socket` fill one socket in order. An explicit `OMP_PLACES` list names the operating system's core numbers, which is why the lists below step by 2.
- **A run smaller than a socket fills socket 0 first.** Runs of 64 cores or more use 64 cores of each node, 32 on each socket, so that every layout of experiment 3 fits the sockets.

| Runs | Launcher options | OpenMP places |
| --- | --- | --- |
| Both serial references (`--launch-serial`) and runs of up to 32 ranks × 1 thread (experiments 1 and 5) | `--map-by core --bind-to core` | `cores` |
| 64 ranks or more (experiments 1, 4 and 5) and the first six layouts of experiment 3 (256 × 1 to 8 × 32) | `--map-by ppr:{ranks_per_socket}:socket:PE={threads} --bind-to core`, with `--cores-per-node 64 --sockets-per-node 2`: 32 cores per socket | `cores` |
| Experiment 2, 1 rank × T threads | the rank on the whole node: `--map-by ppr:1:node --bind-to none` | `"{0}:36:2,{1}:36:2"`: the first 36 threads on socket 0, then socket 1 |
| Experiment 3, 4 × 64, in a runner call of its own (`--run p4_t64_2d_blocking`) | the rank on the whole node: `--map-by ppr:1:node --bind-to none` | `"{0}:32:2,{1}:32:2"`: 32 cores of each socket |

`OMP_PROC_BIND=close` in all cases. Set an explicit `OMP_PLACES` in the environment of that runner call and pass it on with `--mpi-arg=-x --mpi-arg=OMP_PLACES` (the runner reads braces in `--mpi-arg` as placeholders, so the list itself must not go there). Open MPI 5 refuses to bind one rank to cores of both sockets (`PE=72`), so a rank that holds a whole node is not bound by `mpirun`; the explicit `OMP_PLACES` list pins its threads, and the solver's binding check confirms one thread per core. The solvers record the number of nodes; mark in the figures where the runs start to span nodes (above 64 cores).

## GPU build

For experiment 3's GPU runs, build `heat_parallel` with OpenMP offload on a GPU node, in a build folder of its own:

```bash
module load nvhpc-hpcx-cuda12/26.1 CMake/4.0.3-GCCcore-14.3.0
cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=nvc++ -DHEAT2D_OFFLOAD=ON
cmake --build build-gpu -j
```

`HEAT2D_OFFLOAD_FLAGS` defaults to `-mp=gpu;-gpu=cc86,nofma;-Mnofma` for `nvc++`; for `clang++` (`LLVM/20.1.7`) it is `--offload-arch=sm_86`. Untested so far: whether CMake finds HPC-X's MPI and `nvc++`'s OpenMP this way. Before timing, run the blocking cases of experiment 0 with this build (`scripts/check_correctness.py --build-dir build-gpu` with the blocking `--run` names) and check that every run reports `poisson_device` `gpu`. Pass `--build-dir build-gpu` to the runner for the GPU runs.

## Experiment jobs

[jobs/](jobs) holds one PBS script per part of the experiments; [jobs/common.sh](jobs/common.sh) loads the modules, builds in a folder of the job's own, runs the correctness cases as a gate and sets the placement above. Submit from the repository folder, with the queue and project on the command line:

```bash
qsub -q auto_free hpc/jobs/exp0_correctness.pbs
qsub -P <project> hpc/jobs/exp1_one_node.pbs
qsub -P <project> hpc/jobs/exp2_openmp.pbs
qsub -P <project> hpc/jobs/exp5_one_node.pbs
qsub -P <project> hpc/jobs/four_nodes.pbs
qsub -P <project> hpc/jobs/exp3_gpu.pbs
qsub -P <project> hpc/jobs/extra_multinode.pbs
qsub -P <project> hpc/jobs/extra_onenode.pbs
qsub -P <project> hpc/jobs/extra_threads.pbs
qsub -P <project> hpc/jobs/cpu_sizes.pbs
qsub -P <project> hpc/jobs/gpu_sizes.pbs
qsub -P <project> hpc/jobs/exp3d_one_node.pbs
qsub -P <project> hpc/jobs/exp3d_four_nodes.pbs
```

With `-P` and no `-q`, the default queue `auto` routes a job to `batch_cpu` or `batch_gpu`, which accept jobs only from it.

`exp3d_one_node.pbs` and `exp3d_four_nodes.pbs` run the 3D pigment experiments (`experiments/*_3d.json`; `docs/design_3d.md`, model A), each behind the 3D correctness gate (`gate3d` in `common.sh`), with the same placement as their 2D counterparts. `four_nodes.pbs` runs everything above one node (experiments 1 and 5 at 128 and 256 ranks, 3 and 4); `exp3_gpu.pbs` builds the GPU version. The follow-up runs ([experiments/README.md](../experiments/README.md)) are submitted from a checkout of their own, because a job builds what its checkout holds when it starts. `extra_multinode.pbs` (four nodes) runs experiment 5 at 144 ranks, experiment 3 again with 10 repetitions, and 64 ranks on 1, 2 and 4 nodes; `extra_onenode.pbs` runs experiment 5 at 9 and 36 ranks, the placement control with the CPU clock log, and the 64 cores of a node split between ranks and threads; `extra_threads.pbs` runs experiment 2 again; `gpu_sizes.pbs` and `cpu_sizes.pbs` run the GPU and the CPU side of the crossover at four grid sizes. A failed run is recorded and the runner goes on with the next configuration; each run's launcher output, with `--report-bindings`, is in `outputs/<experiment>/<session>/`. Do not change the checkout on Vanda while jobs are queued: each job builds what is there when it starts. Vanda's Open MPI 5 starts its daemons over ssh, even on the job's own node, and they are found only with `mpirun --prefix $EBROOTPRRTE`, and the ranks find the modules' libraries only with `-x LD_LIBRARY_PATH`; [jobs/mpirun.sh](jobs/mpirun.sh) adds both, and the jobs use it as the launcher.

## Large single runs

[big_run.sh](big_run.sh) runs one large case with snapshots, with `heat_serial` or with `heat_parallel` on `RANKS` × `THREADS` cores; its settings come from environment variables. Add the scheduler directives for the platform and set `LAUNCHER` (for example `srun`) before submitting it on a cluster. The main README, section "Larger runs on a bigger machine", explains how the cost grows with the grid. Such runs are for looking at the dynamics; the timed measurements for the report follow the run sequence below.

## Run sequence after implementation

1. Build in Release mode with the selected compiler and MPI environment.
2. Run the correctness checks (experiment 0, `scripts/check_correctness.py`) inside an allocation. The timings are only believed if they pass with this build on this machine.
3. Verify the CPU allocation, actual thread count and rank placement. The run report and the `RESULT` line give the ranks, threads and nodes the solver actually used, the fewest physical cores any rank may run on (`cores_per_rank`) and the `binding`. `heat_parallel` warns when a rank has fewer cores than threads (the launcher has bound the ranks too tightly; with Open MPI, for example, one or two ranks are bound to a single core each by default) and when the binding is not `cores`. Follow "Binding" above until it is `cores`, and record the options used.
4. Run experiments 1 to 5 of [experiments/README.md](../experiments/README.md), in that order, with the planned warmups and measured repetitions and without unrelated concurrent work. Time one run of the benchmark problem first, and adjust its size in `experiments/benchmark_problem.json` if needed. Each experiment runs with [scripts/run_benchmarks.py](../scripts/run_benchmarks.py) inside the job: check the plan with `--dry-run` first, then run it with the launcher options found under "Binding", `--launch-serial` and `--require-binding`, and `--allocated-cores` for a job on several nodes. Placement placeholders such as `{ranks_per_socket}` are described in `scripts/README.md`. Open MPI may not pass `OMP_PLACES` and `OMP_PROC_BIND` on to ranks on other nodes; add `-x OMP_PLACES -x OMP_PROC_BIND` if the first multi-node run reports `binding` `unpinned`.
5. Copy compact raw metrics and environment records into `results/raw/`.
6. Generate figures from those metrics and cite the exact commit/configuration in the report.

Record the scheduler job ID, launcher, module list and hardware details with the results. Large solver snapshots belong in temporary output storage, not Git.
