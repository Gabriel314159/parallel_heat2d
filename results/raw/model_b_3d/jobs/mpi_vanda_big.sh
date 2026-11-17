#!/usr/bin/env bash
# CPU launcher for whole nodes (one or more), one thread per rank, placed as the
# 2D experiments are: up to 32 ranks by core, more ranks 32 per socket (64 per node).
# Use only inside a PBS allocation. MODELB_DRY_RUN=1 prints the command instead.
set -euo pipefail
: "${PBS_NODEFILE:?Run this launcher inside a PBS job}"
threads="${OMP_NUM_THREADS:-1}"
ranks=""
args=("$@")
for ((i = 0; i + 1 < ${#args[@]}; ++i)); do
  if [[ "${args[i]}" == -np ]]; then ranks="${args[i+1]}"; fi
done
if [[ ! "$ranks" =~ ^[1-9][0-9]*$ || ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "mpi_vanda_big.sh: needs -np RANKS and a valid OMP_NUM_THREADS" >&2
  exit 2
fi
allocated=$(sort -u "$PBS_NODEFILE" | wc -l)
if [[ "${MODELB_VALIDATION:-0}" == 1 ]]; then
  # Small correctness runs: all on the first node, unbound.
  placement=(--oversubscribe --map-by slot --bind-to none)
else
  if [[ "$threads" != 1 ]]; then
    echo "mpi_vanda_big.sh: timed runs use one thread per rank" >&2
    exit 2
  fi
  nodes=$(( (ranks + 63) / 64 ))
  if (( nodes > allocated )); then
    echo "mpi_vanda_big.sh: $ranks ranks need $nodes nodes, $allocated allocated" >&2
    exit 2
  fi
  if (( ranks <= 32 )); then
    placement=(--map-by core --bind-to core)
  elif (( ranks % (2 * nodes) == 0 )); then
    placement=(--map-by "ppr:$(( ranks / (2 * nodes) )):socket" --bind-to core)
  else
    echo "mpi_vanda_big.sh: $ranks ranks do not split evenly over $nodes nodes of two sockets" >&2
    exit 2
  fi
fi
prefix=()
if [[ -n "${EBROOTPRRTE:-}" ]]; then prefix=(--prefix "$EBROOTPRRTE"); fi
command=(mpirun "${prefix[@]}" -x LD_LIBRARY_PATH -x OMP_NUM_THREADS -x OMP_DYNAMIC
         -x OMP_PROC_BIND -x OMP_PLACES --report-bindings "${placement[@]}" "$@")
if [[ "${MODELB_DRY_RUN:-0}" == 1 ]]; then
  echo "${command[*]}"
  exit 0
fi
exec "${command[@]}"
