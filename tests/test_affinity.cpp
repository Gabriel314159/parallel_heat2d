// Unit test for src/common/affinity.hpp: counting the cores this process may
// run on, finding the core each thread is bound to, and judging the binding of
// a node from those cores. On Linux it binds itself to one CPU and back;
// elsewhere the counts must be 0 or empty (unknown). Built and run by CTest:
// ctest --test-dir build

#include "common/affinity.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

}  // namespace

int main() {
    // The binding of a node, from the core of each thread (-1: not bound to one core).
    check(node_binding({}) == Binding::unknown, "no cores reported: unknown");
    check(node_binding({3}) == Binding::cores, "one thread on its own core: cores");
    check(node_binding({4, 0, 2, 6}) == Binding::cores, "four threads on four cores: cores");
    check(node_binding({0, 2, 0}) == Binding::shared, "two threads on one core: shared");
    check(node_binding({-1, 5}) == Binding::unpinned, "a thread free to run on several cores: unpinned");
    check(node_binding({-1, -1}) == Binding::unpinned, "two free threads are unpinned, not shared");
    check(node_binding({-1, 2, 2}) == Binding::shared, "sharing outranks a free thread");
    check(binding_name(Binding::cores) == "cores" && binding_name(Binding::shared) == "shared" &&
              binding_name(Binding::unpinned) == "unpinned" && binding_name(Binding::unknown) == "unknown",
          "the names of the bindings");

#if defined(__linux__)
    const int cores = allowed_cores();
    check(cores >= 1, "this process may run on at least one core (" + std::to_string(cores) + ")");
    const std::vector<int> free = thread_cores();
    check(free.size() == 1 && (cores > 1 ? free[0] == -1 : free[0] >= 0),
          "unbound, its one thread is free to run on all its cores");

    // Bound to its first allowed CPU, the process has exactly one core.
    cpu_set_t saved;
    CPU_ZERO(&saved);
    const bool read = sched_getaffinity(0, sizeof(saved), &saved) == 0;
    int first = 0;
    while (read && first < CPU_SETSIZE && !CPU_ISSET(first, &saved)) ++first;
    cpu_set_t one;
    CPU_ZERO(&one);
    CPU_SET(first, &one);
    const bool bound = read && sched_setaffinity(0, sizeof(one), &one) == 0;
    check(bound && allowed_cores() == 1, "bound to one CPU, it has exactly one core");
    const std::vector<int> pinned = thread_cores();
    check(bound && pinned.size() == 1 && pinned[0] == affinity_detail::core_of(first) && pinned[0] <= first,
          "bound to one CPU, its thread is bound to that CPU's core");
    if (read) sched_setaffinity(0, sizeof(saved), &saved);
    check(allowed_cores() == cores, "unbound again, it has all its cores back");
#else
    check(allowed_cores() == 0, "0 (unknown) where the system does not report the binding");
    check(thread_cores().empty(), "no cores (unknown) where the system does not report the binding");
#endif

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
