#!/usr/bin/env bash
# Sourced by the big_*.pbs jobs: modules, build, correctness gate, and one
# function per kind of run. Every run has its own result folder under results/<job id>/.
set -euo pipefail
: "${PBS_O_WORKDIR:?Submit with qsub from the model_b_pressure_3d folder}"
: "${PBS_JOBID:?Missing PBS job ID}"
cd "$PBS_O_WORKDIR"
module purge
module load GCC/14.3.0 OpenMPI/5.0.8-GCC-14.3.0 CMake/4.0.3-GCCcore-14.3.0
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export MODELB_BUILD_JOBS=16
BUILD="$PBS_O_WORKDIR/build-$PBS_JOBID"
RESULTS="$PBS_O_WORKDIR/results/$PBS_JOBID"
mkdir -p "$RESULTS"
echo "job $PBS_JOBID on $(sort -u "$PBS_NODEFILE" | tr '\n' ' ')"
bash ./build.sh "$BUILD" > "$RESULTS/build.log" 2>&1 || { echo "build failed, see $RESULTS/build.log"; exit 1; }
LAUNCHER="bash $PBS_O_WORKDIR/mpi_vanda_big.sh"

# No timings unless this build passes the 69 checks.
gate() {
  (export MODELB_VALIDATION=1 OMP_PROC_BIND=false
   python3 verify.py --build-dir "$BUILD" --output "$RESULTS/validation" --launcher "$LAUNCHER") \
    > "$RESULTS/validation.log" 2>&1 || { echo "correctness gate failed, see $RESULTS/validation.log"; exit 1; }
  tail -1 "$RESULTS/validation.log"
}

# point NAME RANKS N STEPS WARMUPS REPEATS TIMEOUT_S: an N^3 grid on the unit cube,
# RANKS ranks of one thread. A failed run is reported and the job goes on.
point() {
  local name=$1 ranks=$2 n=$3 steps=$4 warmups=$5 repeats=$6 timeout=$7
  echo "$(date +%H:%M:%S) $name: ${n}^3, $ranks ranks, $steps steps"
  if python3 run_model_b.py --suite benchmark --build-dir "$BUILD" --output "$RESULTS/$name" \
       --launcher "$LAUNCHER" --layouts "${ranks}x1" --nx "$n" --ny "$n" --nz "$n" --lx 1 --ly 1 --lz 1 \
       --steps "$steps" --spinup 0 --warmups "$warmups" --repeats "$repeats" --tol 1e-8 \
       --progress 100 --timeout "$timeout" > "$RESULTS/$name.log" 2>&1; then
    tail -1 "$RESULTS/$name.log"
  else
    echo "FAILED: $name"; tail -3 "$RESULTS/$name.log"
  fi
}
