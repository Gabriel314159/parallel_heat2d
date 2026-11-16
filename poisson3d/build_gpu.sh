#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
build="${1:-$root/build-gpu}"
mkdir -p "$build"
flags=(-std=c++17 -O3 -DNDEBUG -Mnofma)
export OMPI_CXX="${NVHPC_CXX:-nvc++}"
"${MPICXX:-mpicxx}" "${flags[@]}" -mp=multicore -DUSE_MPI "$root/poisson3d.cpp" -o "$build/poisson3d_omp"
"${MPICXX:-mpicxx}" "${flags[@]}" -mp=gpu -gpu=cc86,nofma -DUSE_GPU -DUSE_MPI \
    "$root/poisson3d.cpp" -o "$build/poisson3d_gpu"
