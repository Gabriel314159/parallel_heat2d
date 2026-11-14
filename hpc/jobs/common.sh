# Sourced by the CPU jobs: modules, build, gate, placement (hpc/README.md).

cd "$PBS_O_WORKDIR" || exit 1
module purge
module load GCC/14.3.0 OpenMPI/5.0.8-GCC-14.3.0 CMake/4.0.3-GCCcore-14.3.0
export OMP_PLACES=cores OMP_PROC_BIND=close
MPIRUN="$PBS_O_WORKDIR/hpc/jobs/mpirun.sh"
echo "job $PBS_JOBID on $(sort -u "$PBS_NODEFILE" | tr '\n' ' ')"

# A build folder per job: jobs may run at once.
BUILD="build-$PBS_JOBNAME-${PBS_JOBID%%.*}"
cmake -S . -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DHEAT2D_NATIVE=ON -DHEAT2D_BUILD_TESTS="${HEAT2D_TESTS:-OFF}" \
    -DMPIEXEC_EXECUTABLE="$MPIRUN" \
    > "$BUILD.log" 2>&1 && cmake --build "$BUILD" -j 32 >> "$BUILD.log" 2>&1 ||
    { echo "build failed, see $BUILD.log"; exit 1; }

# Correctness gate: no timings unless this build passes.
gate() {
    (unset OMP_PLACES OMP_PROC_BIND
     python3 scripts/check_correctness.py --build-dir "$BUILD" --output-dir "outputs/gate_$PBS_JOBID" \
         --launcher "$MPIRUN --oversubscribe --bind-to none -n {ranks}") > "$BUILD.gate.log" 2>&1 ||
        { echo "correctness gate failed, see $BUILD.gate.log"; exit 1; }
    tail -1 "$BUILD.gate.log"
}

# The same gate for the 3D pigment model (experiments/correctness_3d.json).
gate3d() {
    (unset OMP_PLACES OMP_PROC_BIND
     python3 scripts/check_correctness.py --experiment experiments/correctness_3d.json --build-dir "$BUILD" \
         --output-dir "outputs/gate3d_$PBS_JOBID" \
         --launcher "$MPIRUN --oversubscribe --bind-to none -n {ranks}") > "$BUILD.gate3d.log" 2>&1 ||
        { echo "3D correctness gate failed, see $BUILD.gate3d.log"; exit 1; }
    tail -1 "$BUILD.gate3d.log"
}

RUN=(python3 scripts/run_benchmarks.py --build-dir "$BUILD" --launcher "$MPIRUN" --launch-serial --require-binding
     --allocated-cores "$(wc -l < "$PBS_NODEFILE")" --timeout "${HEAT2D_TIMEOUT:-3600}" --mpi-arg=--report-bindings
     --mpi-arg=-x --mpi-arg=OMP_PLACES --mpi-arg=-x --mpi-arg=OMP_PROC_BIND)
CORE=(--mpi-arg=--map-by --mpi-arg=core --mpi-arg=--bind-to --mpi-arg=core)
SOCKET=(--mpi-arg=--map-by "--mpi-arg=ppr:{ranks_per_socket}:socket:PE={threads}" --mpi-arg=--bind-to --mpi-arg=core
        --cores-per-node 64 --sockets-per-node 2)
# One rank per node: Open MPI binds no rank across sockets, so OpenMP's places pin the threads.
NODE=(--mpi-arg=--map-by --mpi-arg=ppr:1:node --mpi-arg=--bind-to --mpi-arg=none)

# One runner call; on failure, report and go on.
run() { "${RUN[@]}" "$@" || echo "runner failed: $*"; }
