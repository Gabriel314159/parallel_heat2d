cmake_minimum_required(VERSION 3.16)
project(ParallelHeat2D VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimised build unless asked otherwise (Release adds -O3 -DNDEBUG, which
# also switches off the bounds asserts in grid.hpp).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(HEAT2D_ENABLE_MPI "Build the parallel MPI target once its sources exist" ON)
option(HEAT2D_ENABLE_OPENMP "Enable OpenMP in the parallel target" ON)
option(HEAT2D_BUILD_TESTS "Build the unit tests (run them with ctest)" ON)
option(HEAT2D_NATIVE "Optimise for this machine's CPU (-march=native); the binary may not run elsewhere" OFF)

# Shared components: no MPI or OpenMP.
add_library(heat2d_common STATIC
    src/common/config.cpp
    src/common/velocity.cpp
    src/common/physics.cpp
    src/common/stability.cpp
    src/common/flow.cpp
    src/common/io.cpp
)
target_include_directories(heat2d_common PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_definitions(heat2d_common PRIVATE HEAT2D_BUILD_TYPE="${CMAKE_BUILD_TYPE}")
# std::filesystem needs an extra library before GCC 9.1.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
    target_link_libraries(heat2d_common PUBLIC stdc++fs)
endif()

add_executable(heat_serial src/serial/main.cpp src/serial/solver.cpp)
target_link_libraries(heat_serial PRIVATE heat2d_common)

set(HEAT2D_TARGETS heat2d_common heat_serial)

# Parallel solver (src/parallel): configured only when its sources exist.
set(HEAT2D_PARALLEL_SOURCES
    src/parallel/main.cpp
    src/parallel/solver.cpp
    src/parallel/decomposition.cpp
    src/parallel/halo.cpp
)
set(HEAT2D_PARALLEL_READY ON)
foreach(HEAT2D_SOURCE IN LISTS HEAT2D_PARALLEL_SOURCES)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${HEAT2D_SOURCE}")
        set(HEAT2D_PARALLEL_READY OFF)
    endif()
endforeach()

if(HEAT2D_ENABLE_MPI AND NOT HEAT2D_PARALLEL_READY)
    message(STATUS "heat_parallel: sources in src/parallel are not implemented yet; building the serial solver only.")
elseif(HEAT2D_ENABLE_MPI)
    # Optional: without MPI, the serial solver is still built.
    find_package(MPI COMPONENTS CXX)
endif()

if(HEAT2D_ENABLE_MPI AND HEAT2D_PARALLEL_READY AND NOT MPI_CXX_FOUND)
    message(WARNING "MPI not found: heat_parallel is not built. Install an MPI library (e.g. Open MPI) or load its module, then reconfigure.")
elseif(HEAT2D_ENABLE_MPI AND HEAT2D_PARALLEL_READY)
    add_executable(heat_parallel ${HEAT2D_PARALLEL_SOURCES})
    target_link_libraries(heat_parallel PRIVATE heat2d_common MPI::MPI_CXX)

    if(HEAT2D_ENABLE_OPENMP)
        # Apple's clang has no OpenMP runtime of its own; use Homebrew's libomp if present.
        if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND NOT OpenMP_CXX_FLAGS)
            execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE HEAT2D_LIBOMP
                            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
            if(HEAT2D_LIBOMP AND EXISTS "${HEAT2D_LIBOMP}/include/omp.h")
                set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I${HEAT2D_LIBOMP}/include" CACHE STRING "" FORCE)
                set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
                set(OpenMP_omp_LIBRARY "${HEAT2D_LIBOMP}/lib/libomp.dylib" CACHE FILEPATH "" FORCE)
            endif()
        endif()
        find_package(OpenMP COMPONENTS CXX)
        if(OpenMP_CXX_FOUND)
            target_link_libraries(heat_parallel PRIVATE OpenMP::OpenMP_CXX)
            target_compile_definitions(heat_parallel PRIVATE HEAT2D_USE_OPENMP=1)
        else()
            message(WARNING "OpenMP not found: heat_parallel is built MPI-only (one thread per rank).")
        endif()
    endif()
    list(APPEND HEAT2D_TARGETS heat_parallel)
endif()

if(HEAT2D_BUILD_TESTS)
    enable_testing()
    foreach(HEAT2D_TEST velocity config grid io poisson flow)
        add_executable(test_${HEAT2D_TEST} tests/test_${HEAT2D_TEST}.cpp)
        target_link_libraries(test_${HEAT2D_TEST} PRIVATE heat2d_common)
        add_test(NAME ${HEAT2D_TEST} COMMAND test_${HEAT2D_TEST})
        list(APPEND HEAT2D_TARGETS test_${HEAT2D_TEST})
    endforeach()

    add_executable(test_serial_solver tests/test_serial_solver.cpp src/serial/solver.cpp)
    target_link_libraries(test_serial_solver PRIVATE heat2d_common)
    add_test(NAME serial_solver COMMAND test_serial_solver)
    list(APPEND HEAT2D_TARGETS test_serial_solver)

    # End-to-end checks through the command-line interface (scripts/check_correctness.py).
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME correctness_serial
                 COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_correctness.py"
                         --serial-only --build-dir "${CMAKE_BINARY_DIR}"
                         --output-dir "${CMAKE_BINARY_DIR}/correctness"
                 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        # Serial vs parallel comparisons; needs up to 6 MPI ranks. Extra launcher
        # flags (e.g. --oversubscribe on a small machine) go in HEAT2D_MPIEXEC_FLAGS.
        if(TARGET heat_parallel AND MPIEXEC_EXECUTABLE)
            add_executable(test_halo tests/test_halo.cpp src/parallel/decomposition.cpp src/parallel/halo.cpp)
            target_link_libraries(test_halo PRIVATE heat2d_common MPI::MPI_CXX)
            target_compile_options(test_halo PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
            foreach(HEAT2D_RANKS 1 4 6)
                add_test(NAME halo_${HEAT2D_RANKS}_ranks
                         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${HEAT2D_RANKS} ${HEAT2D_MPIEXEC_FLAGS}
                                 $<TARGET_FILE:test_halo>)
            endforeach()
            set(HEAT2D_MPIEXEC_FLAGS "" CACHE STRING "Extra flags for mpiexec in the parallel correctness test")

            # How heat_parallel picks its threads per rank: --threads, else
            # OMP_NUM_THREADS, else 1. Tiny runs; the RESULT line must report
            # the number of threads actually in use.
            set(HEAT2D_THREAD_RUN ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${HEAT2D_MPIEXEC_FLAGS}
                $<TARGET_FILE:heat_parallel> --nx 16 --ny 32 --steps 2 --spinup 0.01)
            if(OpenMP_CXX_FOUND)
                add_test(NAME threads_flag COMMAND ${HEAT2D_THREAD_RUN} --threads 3)
                set_tests_properties(threads_flag PROPERTIES PASS_REGULAR_EXPRESSION "\"omp_threads\": 3,")
                add_test(NAME threads_flag_beats_env COMMAND ${HEAT2D_THREAD_RUN} --threads 2)
                set_tests_properties(threads_flag_beats_env PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4"
                                     PASS_REGULAR_EXPRESSION "\"omp_threads\": 2,")
                add_test(NAME threads_from_env COMMAND ${HEAT2D_THREAD_RUN})
                set_tests_properties(threads_from_env PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=3"
                                     PASS_REGULAR_EXPRESSION "\"omp_threads\": 3,")
                add_test(NAME threads_default COMMAND ${CMAKE_COMMAND} -E env --unset=OMP_NUM_THREADS ${HEAT2D_THREAD_RUN})
                set_tests_properties(threads_default PROPERTIES PASS_REGULAR_EXPRESSION "\"omp_threads\": 1,")
            else()
                add_test(NAME threads_without_openmp COMMAND ${HEAT2D_THREAD_RUN} --threads 2)
                set_tests_properties(threads_without_openmp PROPERTIES PASS_REGULAR_EXPRESSION "built without OpenMP")
            endif()

            add_test(NAME correctness_parallel
                     COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_correctness.py"
                             --build-dir "${CMAKE_BINARY_DIR}"
                             --output-dir "${CMAKE_BINARY_DIR}/correctness"
                             --launcher "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} {ranks} ${HEAT2D_MPIEXEC_FLAGS}"
                     WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        endif()
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(HEAT2D_CPU_FLAGS "")
    if(HEAT2D_NATIVE)
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag(-march=native HEAT2D_HAS_MARCH_NATIVE)
        check_cxx_compiler_flag(-mcpu=native HEAT2D_HAS_MCPU_NATIVE)
        if(HEAT2D_HAS_MARCH_NATIVE)
            set(HEAT2D_CPU_FLAGS -march=native)
        elseif(HEAT2D_HAS_MCPU_NATIVE)
            set(HEAT2D_CPU_FLAGS -mcpu=native)
        else()
            message(WARNING "HEAT2D_NATIVE: the compiler accepts neither -march=native nor -mcpu=native")
        endif()
    endif()
    foreach(HEAT2D_TARGET IN LISTS HEAT2D_TARGETS)
        # -ffp-contract=off: never fuse a*b+c into one rounding step. Every
        # build then rounds each flux the same way, so the serial and parallel
        # solvers can agree bit for bit and both cells sharing a face compute
        # the identical flux.
        target_compile_options(${HEAT2D_TARGET} PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off ${HEAT2D_CPU_FLAGS})
    endforeach()
endif()
