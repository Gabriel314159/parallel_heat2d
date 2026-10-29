// The cores a process may run on, to check how the launcher bound the ranks
// and threads to cores. See docs/design.md, "Binding check".
//
// Header-only, like poisson.hpp: the parallel solver compiles it with OpenMP
// (HEAT2D_USE_OPENMP), so that every thread of the team is asked.

#pragma once

#include <fstream>
#include <set>
#include <string>

#if defined(__linux__)
#include <sched.h>
#endif

// The physical cores this process's threads may run on: the union of the
// affinity masks of all OpenMP threads, with the hardware threads of one core
// counted once. Call it after the thread count is set. Linux only; 0 where the
// system does not tell (e.g. macOS).
inline int allowed_cores() {
#if defined(__linux__)
    cpu_set_t all;
    CPU_ZERO(&all);
    bool known = true;
#ifdef HEAT2D_USE_OPENMP
#pragma omp parallel
#endif
    {
        cpu_set_t mine;
        CPU_ZERO(&mine);
        const bool ok = sched_getaffinity(0, sizeof(mine), &mine) == 0;  // 0: the calling thread
#ifdef HEAT2D_USE_OPENMP
#pragma omp critical(heat2d_allowed_cores)
#endif
        {
            if (ok) {
                CPU_OR(&all, &all, &mine);
            } else {
                known = false;
            }
        }
    }
    if (!known) return 0;

    // The hardware threads of one core share one list of siblings, e.g. "0,64".
    // Where the topology cannot be read, count the CPU as a core of its own.
    std::set<std::string> cores;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        if (!CPU_ISSET(cpu, &all)) continue;
        std::ifstream in("/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/thread_siblings_list");
        std::string siblings;
        cores.insert(std::getline(in, siblings) ? siblings : "cpu " + std::to_string(cpu));
    }
    return static_cast<int>(cores.size());
#else
    return 0;
#endif
}
