// The cores a process may run on, to check how the launcher bound the ranks
// and threads to cores. See docs/design.md, "Binding check".
//
// Header-only, like poisson.hpp: the parallel solver compiles it with OpenMP
// (HEAT2D_USE_OPENMP), so that every thread of the team is asked.

#pragma once

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif
#ifdef HEAT2D_USE_OPENMP
#include <omp.h>
#endif

// How the threads of a run are bound to physical cores, from best to worst.
enum class Binding {
    cores,     // every thread is bound to a physical core of its own
    unpinned,  // some thread may run on more than one core
    shared,    // two threads are bound to the same core, e.g. its two hardware threads
    unknown    // the system does not tell (e.g. macOS)
};

inline std::string binding_name(Binding b) {
    switch (b) {
        case Binding::cores:
            return "cores";
        case Binding::unpinned:
            return "unpinned";
        case Binding::shared:
            return "shared";
        default:
            return "unknown";
    }
}

// The binding of the threads of one node, given the core each thread is bound
// to (thread_cores; -1 for a thread that may run on more than one core).
// Unknown if the list is empty.
inline Binding node_binding(std::vector<int> cores) {
    if (cores.empty()) return Binding::unknown;
    std::sort(cores.begin(), cores.end());  // the -1 come first
    for (std::size_t k = 1; k < cores.size(); ++k) {
        if (cores[k] >= 0 && cores[k] == cores[k - 1]) return Binding::shared;
    }
    return cores.front() < 0 ? Binding::unpinned : Binding::cores;
}

#if defined(__linux__)
namespace affinity_detail {

// The physical core of a CPU (a hardware thread), named by the lowest CPU
// number among the hardware threads of that core: the first number of its
// list of siblings, e.g. "0,64" or "0-1". Where the topology cannot be read,
// the CPU counts as a core of its own.
inline int core_of(int cpu) {
    std::ifstream in("/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/thread_siblings_list");
    std::string siblings;
    if (!std::getline(in, siblings) || siblings.empty() || !std::isdigit(static_cast<unsigned char>(siblings[0]))) {
        return cpu;
    }
    return std::stoi(siblings);
}

// The physical cores of the CPUs in a mask.
inline std::set<int> cores_in(const cpu_set_t& mask) {
    std::set<int> cores;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (CPU_ISSET(cpu, &mask)) cores.insert(core_of(cpu));
    }
    return cores;
}

// The affinity mask of each OpenMP thread of this process (each thread asks
// for its own), in thread order. Empty if a mask cannot be read. Call it after
// the thread count is set.
inline std::vector<cpu_set_t> thread_masks() {
    std::vector<cpu_set_t> masks(1);
    std::vector<char> read(1, 0);
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel
#endif
    {
#ifdef HEAT2D_USE_OPENMP
#pragma omp single
        {
            masks.resize(static_cast<std::size_t>(omp_get_num_threads()));
            read.assign(masks.size(), 0);
        }
        const std::size_t t = static_cast<std::size_t>(omp_get_thread_num());
#else
        const std::size_t t = 0;
#endif
        CPU_ZERO(&masks[t]);
        read[t] = sched_getaffinity(0, sizeof(cpu_set_t), &masks[t]) == 0;  // 0: the calling thread
    }
    if (std::find(read.begin(), read.end(), 0) != read.end()) return {};
    return masks;
}

}  // namespace affinity_detail
#endif

// The physical cores this process's threads may run on: the union of the
// affinity masks of all OpenMP threads, with the hardware threads of one core
// counted once. Call it after the thread count is set. Linux only; 0 where the
// system does not tell (e.g. macOS).
inline int allowed_cores() {
#if defined(__linux__)
    const std::vector<cpu_set_t> masks = affinity_detail::thread_masks();
    if (masks.empty()) return 0;
    cpu_set_t all;
    CPU_ZERO(&all);
    for (const cpu_set_t& mask : masks) CPU_OR(&all, &all, &mask);
    return static_cast<int>(affinity_detail::cores_in(all).size());
#else
    return 0;
#endif
}

// The physical core each OpenMP thread of this process is bound to, in thread
// order, named as in core_of; -1 for a thread that may run on more than one
// core. Call it after the thread count is set. Linux only; empty where the
// system does not tell (e.g. macOS).
inline std::vector<int> thread_cores() {
    std::vector<int> cores;
#if defined(__linux__)
    for (const cpu_set_t& mask : affinity_detail::thread_masks()) {
        const std::set<int> mine = affinity_detail::cores_in(mask);
        cores.push_back(mine.size() == 1 ? *mine.begin() : -1);
    }
#endif
    return cores;
}
