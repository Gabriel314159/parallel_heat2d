#!/usr/bin/env bash
# Template only. Add the actual scheduler directives and module setup here.
# Run this inside an allocated HPC job, after building the solver.
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
    echo "The parallel solver has not been built at ./build/heat_parallel." >&2
    exit 2
fi

export OMP_NUM_THREADS
mkdir -p logs

"${MPI_LAUNCH_ARGS[@]}" ./build/heat_parallel \
    --nx 1024 --ny 1024 --steps 1000 \
    --case heated_top --cfl 0.8 \
    --decomp 2d --exchange blocking
