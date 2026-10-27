#!/usr/bin/env bash
# One large run with snapshots, on a workstation or a cluster node.
#
# Uses heat_parallel with RANKS MPI ranks x THREADS OpenMP threads, or
# heat_serial when RANKS=1 and THREADS=1. All settings come from the
# environment; the defaults are 2048 x 2048 cells for two flow periods.
#
#   HEAT2D_REPO=$PWD RANKS=16 bash hpc/big_run.sh
#   HEAT2D_REPO=$PWD RANKS=4 THREADS=4 NX=4096 KAPPA=0.000125 STEPS=338400 \
#       INTERVAL=3384 COARSEN=8 bash hpc/big_run.sh
#
# On a cluster, add the scheduler directives for the platform here and set
# LAUNCHER to its approved launcher. For Slurm they look like the lines below;
# these are examples, not verified settings for Atlas or Vanda:
#   #SBATCH --job-name=heat2d
#   #SBATCH --ntasks=16
#   #SBATCH --cpus-per-task=1
#   #SBATCH --time=01:00:00
# with LAUNCHER="srun".
#
# See "Larger runs on a bigger machine" in README.md for grid sizes, run times
# and snapshot sizes.
set -euo pipefail

: "${HEAT2D_REPO:?Set HEAT2D_REPO to the repository path, e.g. HEAT2D_REPO=\$PWD}"
cd "$HEAT2D_REPO"

NX="${NX:-2048}"
KAPPA="${KAPPA:-0.00025}"
STEPS="${STEPS:-169200}"
INTERVAL="${INTERVAL:-1692}"   # steps between snapshots
COARSEN="${COARSEN:-4}"        # snapshots as COARSEN x COARSEN block averages
RANKS="${RANKS:-1}"
THREADS="${THREADS:-1}"
DECOMP="${DECOMP:-2d}"
EXCHANGE="${EXCHANGE:-overlap}"
LAUNCHER="${LAUNCHER:-mpiexec -n}"
RUN_DIR="${RUN_DIR:-outputs/run${NX}}"

PROBLEM=(--nx "$NX" --ny "$NX" --kappa "$KAPPA" --steps "$STEPS"
         --snapshot-interval "$INTERVAL" --snapshot-coarsen "$COARSEN" --snapshot-dir "$RUN_DIR"
         --output "$RUN_DIR/final.npy")

if [[ "$RANKS" == 1 && "$THREADS" == 1 ]]; then
    COMMAND=(./build/heat_serial "${PROBLEM[@]}")
else
    read -r -a LAUNCH <<< "$LAUNCHER"
    COMMAND=("${LAUNCH[@]}" "$RANKS" ./build/heat_parallel "${PROBLEM[@]}" --decomp "$DECOMP" --exchange "$EXCHANGE")
fi

if [[ ! -x "${COMMAND[0]}" ]] && ! command -v "${COMMAND[0]}" > /dev/null; then
    echo "${COMMAND[0]} not found. Build first: cmake -S . -B build && cmake --build build --parallel" >&2
    exit 2
fi
if [[ "$RANKS" != 1 || "$THREADS" != 1 ]] && [[ ! -x ./build/heat_parallel ]]; then
    echo "./build/heat_parallel not found: install MPI (and OpenMP), then rebuild. Or use RANKS=1 THREADS=1." >&2
    exit 2
fi

export OMP_NUM_THREADS="$THREADS"
mkdir -p "$RUN_DIR"
echo "${NX} x ${NX}, kappa ${KAPPA}, ${STEPS} steps, ${RANKS} ranks x ${THREADS} threads -> ${RUN_DIR}"
echo "${COMMAND[*]}"
"${COMMAND[@]}" | tee "$RUN_DIR/report.txt"
