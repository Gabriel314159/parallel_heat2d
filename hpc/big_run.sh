#!/usr/bin/env bash
# One large run with snapshots, on a workstation or a cluster node.
#
# Uses heat_parallel with RANKS MPI ranks x THREADS OpenMP threads, or
# heat_serial when RANKS=1 and THREADS=1. All settings come from the
# environment. The defaults: the 1 x 2 basin on 512 x 1024 cells with the
# default flow (viscosity 1e-3, wind 0.1), a spin-up of 100, then about 40
# time units with the pigment (dt is about 2e-4), 100 snapshots.
#
#   HEAT2D_REPO=$PWD RANKS=16 bash hpc/big_run.sh
#   HEAT2D_REPO=$PWD RANKS=16 NX=1024 VISCOSITY=0.00025 WIND=0.05 MAXSPEED=3 \
#       SPINUP=300 STEPS=600000 INTERVAL=6000 COARSEN=4 bash hpc/big_run.sh
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

NX="${NX:-512}"                 # cells along x; the basin has 2 * NX cells along y
VISCOSITY="${VISCOSITY:-0.001}"
WIND="${WIND:-0.1}"             # wind amplitude F0
MAXSPEED="${MAXSPEED:-2}"       # speed estimate that sets dt
SPINUP="${SPINUP:-100}"
KAPPA="${KAPPA:-0.0001}"
STEPS="${STEPS:-200000}"
INTERVAL="${INTERVAL:-2000}"    # steps between snapshots
COARSEN="${COARSEN:-2}"         # snapshots shrunk COARSEN times in each direction
EXTRA="${EXTRA:-}"              # any further solver flags, e.g. "--wind periodic --wind-alpha 1.5 --wind-period 20"
RANKS="${RANKS:-1}"
THREADS="${THREADS:-1}"
DECOMP="${DECOMP:-2d}"
EXCHANGE="${EXCHANGE:-overlap}"
LAUNCHER="${LAUNCHER:-mpiexec -n}"
RUN_DIR="${RUN_DIR:-outputs/run${NX}}"

read -r -a EXTRA_ARGS <<< "$EXTRA"
PROBLEM=(--nx "$NX" --ny "$((2 * NX))" --viscosity "$VISCOSITY" --wind-amplitude "$WIND" --max-speed "$MAXSPEED"
         --spinup "$SPINUP" --kappa "$KAPPA" --steps "$STEPS"
         --snapshot-interval "$INTERVAL" --snapshot-coarsen "$COARSEN" --snapshot-dir "$RUN_DIR"
         --output "$RUN_DIR/final.npy" "${EXTRA_ARGS[@]}")

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
echo "${NX} x $((2 * NX)), viscosity ${VISCOSITY}, wind ${WIND}, ${STEPS} steps, ${RANKS} ranks x ${THREADS} threads -> ${RUN_DIR}"
echo "${COMMAND[*]}"
"${COMMAND[@]}" | tee "$RUN_DIR/report.txt"
