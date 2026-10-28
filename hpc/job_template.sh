#!/usr/bin/env bash
# Template only. Add the actual scheduler directives and module setup here.
# Run this inside an allocated HPC job, after implementing and building the solver.
# MPI_LAUNCH_ARGS must contain the platform-approved launcher and rank options.
set -euo pipefail

: "${HEAT2D_REPO:?Set HEAT2D_REPO to the repository path visible from the compute node.}"
: "${OMP_NUM_THREADS:?Set OMP_NUM_THREADS to the allocated threads per rank.}"

cd "$HEAT2D_REPO"

# Replace the empty array with the verified launcher. Example syntax varies by platform.
MPI_LAUNCH_ARGS=()
if [[ ${#MPI_LAUNCH_ARGS[@]} -eq 0 ]]; then
    echo "Configure MPI_LAUNCH_ARGS and scheduler resources before using this template." >&2
    exit 2
fi

if [[ ! -x ./build/heat_parallel ]]; then
    echo "The parallel solver has not been implemented/built at ./build/heat_parallel." >&2
    exit 2
fi

export OMP_NUM_THREADS
mkdir -p logs

# Example values only: choose the grid and steps from experiments/*.json.
# --spinup 0 starts the flow at rest, so every run does the same work.
"${MPI_LAUNCH_ARGS[@]}" ./build/heat_parallel \
    --nx 512 --ny 1024 --steps 200 --spinup 0 \
    --kappa 0.0001 --cfl 0.8 --max-speed 2 \
    --viscosity 0.001 --wind steady --wind-amplitude 0.1 \
    --case pigment_drop --blob-x 0.2 --blob-y 1.0 --blob-width 0.05 \
    --decomp 2d --exchange blocking
