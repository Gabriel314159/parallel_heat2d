#!/usr/bin/env bash
set -euo pipefail
: "${PBS_O_WORKDIR:?Submit with qsub from the model_b folder}"
: "${PBS_JOBID:?Missing PBS job ID}"
cd "$PBS_O_WORKDIR"
module purge
module load GCC/14.3.0 OpenMPI/5.0.8-GCC-14.3.0 CMake/4.0.3-GCCcore-14.3.0
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export MODELB_BUILD_JOBS=4
BUILD="$PBS_O_WORKDIR/build-$PBS_JOBID"
RESULTS="$PBS_O_WORKDIR/results/$PBS_JOBID"
mkdir -p "$RESULTS"
trap 'status=$?; if ((status != 0)); then echo "Model B job FAILED (exit $status). Check results/$PBS_JOBID and the PBS log."; fi' EXIT
bash ./build.sh "$BUILD" > "$RESULTS/build.log" 2>&1
LAUNCHER="bash $PBS_O_WORKDIR/mpi_vanda.sh"
