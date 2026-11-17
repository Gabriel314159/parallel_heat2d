#!/usr/bin/env bash
# CPU launcher. Use only inside a PBS allocation. Intended for the supplied
# single-node jobs; a separate allocation/placement audit is needed for multi-node work.
set -euo pipefail
: "${PBS_NODEFILE:?Run this launcher inside a PBS job}"
mapfile -t hosts < <(sort -u "$PBS_NODEFILE")
if [[ ${#hosts[@]} != 1 ]]; then
  echo "This checked launcher expects a single allocated node." >&2
  exit 2
fi
threads="${OMP_NUM_THREADS:-1}"
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid OMP_NUM_THREADS" >&2
  exit 2
fi
prefix=()
if [[ -n "${EBROOTPRRTE:-}" ]]; then
  prefix=(--prefix "$EBROOTPRRTE")
fi
if [[ "${MODELB_VALIDATION:-0}" == 1 ]]; then
  # The correctness gate can use 8 ranks even if PBS mpiprocs is 4 for a later
  # 4x8 hybrid run. All supplied allocations have at least 16 real CPU cores.
  placement=(--oversubscribe --map-by slot --bind-to none)
else
  placement=(--map-by "slot:PE=$threads" --bind-to core --report-bindings)
fi
exec mpirun "${prefix[@]}" -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -x OMP_DYNAMIC \
  -x OMP_PROC_BIND -x OMP_PLACES "${placement[@]}" "$@"
