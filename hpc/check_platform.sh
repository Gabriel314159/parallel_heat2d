#!/usr/bin/env bash
# Report what a cluster offers for our runs: the CPUs, NUMA layout and GPUs of
# the node it runs on, the scheduler's queues and node types, the modules
# (compilers, MPI, CUDA, NVHPC), and which compilers build OpenMP offload code.
# Read-only: it submits no job, loads modules only inside this script, and
# compiles in a temporary directory that it removes at the end.
#
# On a login node, from the laptop (one login; the report is saved locally):
#   ssh <user>@vanda.nus.edu.sg 'bash -l -s' < hpc/check_platform.sh > outputs/vanda_platform.txt
# Vanda's login node is a different machine from its compute nodes, so the
# compile tests run only inside a job. There the script reports the compute
# node, and the offload test shows whether the target region really runs on
# that node's GPU:
#   bash hpc/check_platform.sh > platform_node.txt

progress() { printf 'check_platform: %s\n' "$*" >&2; }
section() { printf '\n===== %s =====\n' "$*"; }
have() { command -v "$1" > /dev/null 2>&1; }

# A time limit where `timeout` exists (not on macOS).
limit() {
    local seconds=$1
    shift
    if have timeout; then timeout "$seconds" "$@"; else "$@"; fi
}

# Print a command line, then run it in this shell (so that `module` works),
# with stderr folded into stdout. Stdin is closed: when the script itself comes
# through stdin (bash -s), a command reading stdin would swallow the rest.
show() {
    printf '$ %s\n' "$1"
    eval "$1" < /dev/null 2>&1
    local status=$?
    if [ "$status" -ne 0 ]; then printf '(exit status %s)\n' "$status"; fi
    return 0
}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
shopt -s expand_aliases # site helpers such as `hpc` may be aliases
in_job=no
if [ -n "${PBS_JOBID:-}" ] || [ -n "${SLURM_JOB_ID:-}" ]; then in_job=yes; fi

# A non-interactive shell may not define `module`.
if ! type module > /dev/null 2>&1; then
    for init in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh \
        /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do
        if [ -r "$init" ]; then
            . "$init" > /dev/null 2>&1
            break
        fi
    done
fi
have_module() { type module > /dev/null 2>&1; }

echo "Platform report by hpc/check_platform.sh"

progress "host, CPUs and GPUs"
section "Host"
show 'date -u "+%Y-%m-%d %H:%M UTC"'
show 'hostname -f 2> /dev/null || hostname'
show 'uname -srm'
if [ -r /etc/os-release ]; then show 'grep -E "^PRETTY_NAME=" /etc/os-release'; fi
show 'env | grep -E "^(PBS_JOBID|PBS_QUEUE|SLURM_JOB_ID|SLURM_JOB_PARTITION|SLURM_CPUS_ON_NODE|OMP_[A-Z_]*)=" || echo "(not inside a job; no OMP_ variables set)"'

section "Site help (the hpc command: sample job scripts, projects, storage)"
if type hpc > /dev/null 2>&1; then
    show 'hpc cpu'
    show 'hpc gpu'
    show 'hpc project'
    show 'hpc s'
else
    echo "(no hpc command)"
fi

section "CPUs and memory of this node"
if have lscpu; then show 'lscpu'; else echo "(no lscpu)"; fi
if have nproc; then show 'nproc'; fi
if [ -r /sys/devices/system/cpu/smt/active ]; then
    printf 'SMT (hyper-threading) active: %s (1 = on)\n' "$(cat /sys/devices/system/cpu/smt/active)"
fi
if have numactl; then show 'numactl --hardware'; fi
if have free; then show 'free -g'; fi

section "GPUs of this node"
if have nvidia-smi; then
    show 'nvidia-smi -L'
    show 'nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv'
    show 'nvidia-smi topo -m'
else
    echo "(no nvidia-smi on this node)"
fi
if have lspci; then
    show 'lspci | grep -iE "nvidia|3d controller|amd/ati" || echo "(no GPU on the PCI bus)"'
fi

progress "scheduler, queues and node types"
section "Scheduler"
scheduler=none
if have qstat && have pbsnodes; then scheduler=pbs; fi
if have sinfo; then scheduler=slurm; fi
echo "scheduler: $scheduler"
case "$scheduler" in
pbs)
    show 'qstat --version'
    show 'qstat -Q'
    show 'qstat -Qf | grep -E "^Queue:|queue_type|enabled|started|resources_(max|min|default|available)|default_chunk|max_run|max_queued|acl_user_enable|route_destinations"'
    show 'limit 60 pbsnodes -aSj'
    echo
    echo "Node types (count, then each node's resources; pcpus counts hardware threads):"
    limit 60 pbsnodes -a < /dev/null 2>&1 | awk '
        function flush() {
            if (name != "") print "ncpus=" c ", pcpus=" pc ", mem=" m ", ngpus=" g ", Qlist=" q other
            name = ""; c = pc = m = g = q = other = ""
        }
        /^[^[:space:]]/ { flush(); name = $1; next }
        $1 == "pcpus" { pc = $3 }
        $1 == "resources_available.ncpus" { c = $3 }
        $1 == "resources_available.mem" { m = $3 }
        $1 == "resources_available.ngpus" { g = $3 }
        tolower($1) == "resources_available.qlist" { q = $3 }
        # site resources that describe the hardware, such as a GPU or CPU model
        $1 ~ /^resources_available\.(arch|model|.*(cpu|gpu|type|feature|socket|numa).*)$/ &&
            $1 != "resources_available.ncpus" && $1 != "resources_available.ngpus" {
            other = other ", " substr($1, 21) "=" $3
        }
        END { flush() }' | sort | uniq -c | sort -rn
    ;;
slurm)
    show 'sinfo --version'
    show 'sinfo -s'
    show "sinfo -o '%P %a %l %D %c %m %G %f'"
    echo
    echo "Node types (count; CPUs, sockets, cores per socket, threads per core, memory in MB, GRES, partition):"
    sinfo -N -h -o '%c %X %Y %Z %m %G %P' < /dev/null 2>&1 | sort | uniq -c | sort -rn
    show 'scontrol show partition | grep -E "PartitionName|MaxNodes|MaxTime|DefaultTime|MaxCPUsPerNode|OverSubscribe|DefMemPerCPU|TRES="'
    show 'scontrol show config | grep -iE "SLURM_VERSION|SelectType|TaskPlugin|TaskAffinity|MpiDefault|GresTypes"'
    show 'sacctmgr -n -P show assoc user="$USER" format=account,partition,qos'
    ;;
*)
    echo "(neither PBS nor Slurm commands are on the PATH)"
    ;;
esac

progress "modules"
section "Modules"
touch "$work/modules.txt"
if have_module; then
    show 'module --version'
    show 'module list'
    module -t avail < /dev/null 2>&1 | grep -v ':$' | sort -u > "$work/modules.txt"
    echo
    echo "Of interest (compilers, MPI, CUDA, NVHPC, CMake, Python):"
    grep -iE 'gcc|gnu|intel|oneapi|openmpi|mpich|mvapich|impi|hpcx|cuda|nvhpc|hpc[-_]sdk|pgi|llvm|clang|aocc|cmake|python|conda' \
        "$work/modules.txt" || echo "(none)"
    if module --version 2>&1 | grep -qi lmod; then
        echo
        echo "Lmod spider (modules that appear only after a compiler or MPI is loaded):"
        module -t spider < /dev/null 2>&1 | grep -iE 'nvhpc|cuda|openmpi|mpich|mvapich|impi|hpcx|llvm|clang|gcc|intel|oneapi' | sort -u
    fi
    echo
    echo "All available modules:"
    cat "$work/modules.txt"
else
    echo "(no module command)"
fi

progress "compilers and MPI"
section "Compilers and MPI on the default PATH"
for tool in gcc g++ clang++ icpx nvc++ nvcc mpicc mpicxx mpirun mpiexec srun cmake python3; do
    if path=$(command -v "$tool"); then printf '%-8s %s\n' "$tool" "$path"; fi
done
if have g++; then
    show 'g++ --version | head -n 1'
    show 'g++ -v 2>&1 | grep -E "^(Target|OFFLOAD_TARGET_NAMES)" || echo "(this g++ was built without offload targets)"'
fi
if have mpirun; then show 'mpirun --version 2>&1 | head -n 2'; fi
if have ompi_info; then
    show 'ompi_info 2> /dev/null | grep -E "Open MPI:|Thread support"'
    show 'ompi_info --parsable --all 2> /dev/null | grep -E "mpi_built_with_cuda_support:value"'
fi
if have cmake; then show 'cmake --version | head -n 1'; fi

progress "OpenMP offload test (a few small compiles)"
section "OpenMP offload"
cat > "$work/offload.cpp" << 'EOF'
#include <cstdio>
#include <omp.h>
int main() {
    const int n = 1 << 20;
    double sum = 0.0;
    int initial = 1;
#pragma omp target teams distribute parallel for reduction(+ : sum) map(tofrom : sum)
    for (int i = 0; i < n; ++i) sum += 1.0;
#pragma omp target map(from : initial)
    { initial = omp_is_initial_device(); }
    std::printf("devices %d; sum %.0f of %d; the target region ran on the %s\n", omp_get_num_devices(), sum, n,
                initial ? "host (no offload)" : "GPU");
    return sum == n ? 0 : 1;
}
EOF
touch "$work/offload_summary.txt"

# try_offload LABEL MODULE COMMAND...: load MODULE ("-" for none), build the
# test with COMMAND, run it, and add one line to the summary.
try_offload() {
    local label=$1 mod=$2
    shift 2
    printf '\n--- %s ---\n' "$label"
    (
        cd "$work" || exit 1
        if [ "$mod" != "-" ] && ! module load "$mod" < /dev/null > load.log 2>&1; then
            cat load.log
            echo "$label: module does not load" >> offload_summary.txt
            exit 1
        fi
        if ! have "$1"; then
            echo "$1 is not on the PATH"
            exit 1
        fi
        printf '$ %s -o offload offload.cpp\n' "$*"
        if ! limit 180 "$@" -o offload offload.cpp < /dev/null > build.log 2>&1; then
            tail -n 12 build.log
            echo "$label: does not build" >> offload_summary.txt
            exit 1
        fi
        tail -n 4 build.log
        limit 60 ./offload < /dev/null > run.log 2>&1
        tail -n 6 run.log
        echo "$label: builds; run: $(tail -n 1 run.log)" >> offload_summary.txt
        rm -f offload
    )
}

newest() { grep -iE "$1" "$work/modules.txt" | sort -V | tail -n 1; }
nvhpc_module=$(newest 'nvhpc|hpc[-_]sdk')
llvm_module=$(newest '^(llvm|clang)')
gcc_module=$(newest '^gcc(/|-)')

if [ "$in_job" = yes ]; then
    try_offload "nvc++ on the PATH" - nvc++ -mp=gpu
    if [ -n "$nvhpc_module" ]; then try_offload "nvc++ from $nvhpc_module" "$nvhpc_module" nvc++ -mp=gpu; fi
    try_offload "clang++ on the PATH" - clang++ -fopenmp --offload-arch=sm_80
    if [ -n "$llvm_module" ]; then
        try_offload "clang++ from $llvm_module" "$llvm_module" clang++ -fopenmp --offload-arch=sm_80
    fi
    try_offload "g++ on the PATH" - g++ -fopenmp -foffload=nvptx-none
    if [ -n "$gcc_module" ]; then
        try_offload "g++ from $gcc_module" "$gcc_module" g++ -fopenmp -foffload=nvptx-none
    fi
else
    echo "(skipped outside a job: compile on a compute node; candidates: ${nvhpc_module:-no NVHPC module}," \
        "${llvm_module:-no LLVM module}, ${gcc_module:-no GCC module})"
    echo "not tested (login node)" > "$work/offload_summary.txt"
fi

progress "storage"
section "Storage"
show 'df -h "$HOME" | tail -n 1'
for dir in /scratch /scratch2 /hpctmp /hpctmp2 "${SCRATCH:-}" "${TMPDIR:-}"; do
    if [ -n "$dir" ] && [ -d "$dir" ]; then show "df -h '$dir' | tail -n 1"; fi
done

section "Summary"
echo "scheduler: $scheduler"
if have lscpu; then
    lscpu | grep -E '^(Model name|Socket\(s\)|Core\(s\) per socket|Thread\(s\) per core|NUMA node\(s\)):'
fi
echo "GPUs on this node: $(nvidia-smi -L 2> /dev/null | grep -c '^GPU')"
echo "OpenMP offload compilers:"
if [ -s "$work/offload_summary.txt" ]; then sed 's/^/  /' "$work/offload_summary.txt"; else echo "  none found"; fi
progress "done"
