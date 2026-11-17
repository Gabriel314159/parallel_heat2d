#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
build="${1:-$root/build}"
args=(-DCMAKE_BUILD_TYPE=Release)
if [[ "${MODELB_SERIAL_ONLY:-0}" == 1 ]]; then
  args+=(-DMODELB_MPI=OFF)
elif [[ "$(uname -s)" == Darwin ]] && command -v brew >/dev/null 2>&1; then
  omp_prefix="$(brew --prefix libomp)"
  args+=("-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I$omp_prefix/include"
         -DOpenMP_CXX_LIB_NAMES=omp "-DOpenMP_omp_LIBRARY=$omp_prefix/lib/libomp.dylib")
fi
cmake -S "$root" -B "$build" "${args[@]}"
cmake --build "$build" -j "${MODELB_BUILD_JOBS:-4}"
