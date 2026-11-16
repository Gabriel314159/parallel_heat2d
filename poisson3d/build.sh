#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
build="${1:-$root/build}"
mkdir -p "$build"
flags=(-std=c++17 -O3 -DNDEBUG -ffp-contract=off)
"${CXX:-g++}" "${flags[@]}" "$root/poisson3d.cpp" -o "$build/poisson3d_serial"
if [[ "${POISSON3D_SERIAL_ONLY:-0}" != 1 ]]; then
    "${CXX:-g++}" "${flags[@]}" -fopenmp "$root/poisson3d.cpp" -o "$build/poisson3d_serial_omp"
    "${MPICXX:-mpicxx}" "${flags[@]}" -DUSE_MPI -fopenmp "$root/poisson3d.cpp" -o "$build/poisson3d_mpi"
fi
