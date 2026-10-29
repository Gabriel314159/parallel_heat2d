// Unit test for src/common/affinity.hpp: counting the cores this process may
// run on. On Linux it binds itself to one CPU and back; elsewhere the count
// must be 0 (unknown). Built and run by CTest: ctest --test-dir build

#include "common/affinity.hpp"

#include <iostream>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& label) {
    std::cout << (condition ? "  ok    " : "  FAIL  ") << label << "\n";
    if (!condition) ++failures;
}

}  // namespace

int main() {
#if defined(__linux__)
    const int cores = allowed_cores();
    check(cores >= 1, "this process may run on at least one core (" + std::to_string(cores) + ")");

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
    if (read) sched_setaffinity(0, sizeof(saved), &saved);
    check(allowed_cores() == cores, "unbound again, it has all its cores back");
#else
    check(allowed_cores() == 0, "0 (unknown) where the system does not report the binding");
#endif

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
