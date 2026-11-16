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
option(HEAT2D_OFFLOAD "heat_parallel: run the CG Poisson solve on a GPU with OpenMP target offload (nvc++ or clang++)" OFF)

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

# The 3D pigment model (docs/design_3d.md, "Model A"): the 3D pigment in the
# 2D flow. Its shared parts add to heat2d_common; no MPI or OpenMP.
add_library(heat3d_common STATIC
    src/common3d/config3d.cpp
    src/common3d/physics3d.cpp
    src/common3d/io3d.cpp
)
target_link_libraries(heat3d_common PUBLIC heat2d_common)

add_executable(heat3d_serial src/serial3d/main.cpp src/serial3d/solver.cpp)
target_link_libraries(heat3d_serial PRIVATE heat3d_common)

set(HEAT2D_TARGETS heat2d_common heat_serial heat3d_common heat3d_serial)

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
endif()

# heat_serial_omp: the serial solver with the OpenMP pragmas, on one thread.
if(HEAT2D_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    add_executable(heat_serial_omp src/serial/main.cpp src/serial/solver.cpp)
    target_link_libraries(heat_serial_omp PRIVATE heat2d_common OpenMP::OpenMP_CXX)
    target_compile_definitions(heat_serial_omp PRIVATE HEAT2D_USE_OPENMP=1 HEAT2D_SERIAL_OMP=1)
    list(APPEND HEAT2D_TARGETS heat_serial_omp)
endif()

# Standalone 3D Poisson benchmark (poisson3d/README.md), built as poisson3d/build.sh does.
add_executable(poisson3d_serial poisson3d/poisson3d.cpp)
list(APPEND HEAT2D_TARGETS poisson3d_serial)
if(OpenMP_CXX_FOUND)
    add_executable(poisson3d_serial_omp poisson3d/poisson3d.cpp)
    target_link_libraries(poisson3d_serial_omp PRIVATE OpenMP::OpenMP_CXX)
    list(APPEND HEAT2D_TARGETS poisson3d_serial_omp)
endif()

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
        if(OpenMP_CXX_FOUND)
            target_link_libraries(heat_parallel PRIVATE OpenMP::OpenMP_CXX)
            target_compile_definitions(heat_parallel PRIVATE HEAT2D_USE_OPENMP=1)
        else()
            message(WARNING "OpenMP not found: heat_parallel is built MPI-only (one thread per rank).")
        endif()
    endif()
    list(APPEND HEAT2D_TARGETS heat_parallel)

    # heat3d_parallel: the 3D pigment model over MPI + OpenMP (slabs or pencils).
    add_executable(heat3d_parallel src/parallel3d/main.cpp src/parallel3d/solver.cpp src/parallel3d/halo3d.cpp
                   src/parallel/decomposition.cpp src/parallel/halo.cpp)
    target_link_libraries(heat3d_parallel PRIVATE heat3d_common MPI::MPI_CXX)
    if(HEAT2D_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(heat3d_parallel PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(heat3d_parallel PRIVATE HEAT2D_USE_OPENMP=1)
    endif()
    list(APPEND HEAT2D_TARGETS heat3d_parallel)

    # poisson3d_mpi: the standalone 3D Poisson benchmark over MPI + OpenMP.
    add_executable(poisson3d_mpi poisson3d/poisson3d.cpp)
    target_link_libraries(poisson3d_mpi PRIVATE MPI::MPI_CXX)
    target_compile_definitions(poisson3d_mpi PRIVATE USE_MPI)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(poisson3d_mpi PRIVATE OpenMP::OpenMP_CXX)
    endif()
    list(APPEND HEAT2D_TARGETS poisson3d_mpi)
endif()

# GPU build: the CG of heat_parallel on an A40 (docs/design.md). Empty
# HEAT2D_OFFLOAD_FLAGS runs the target regions on the host, for testing.
if(HEAT2D_OFFLOAD)
    if(NOT TARGET heat_parallel OR NOT OpenMP_CXX_FOUND)
        message(FATAL_ERROR "HEAT2D_OFFLOAD needs heat_parallel, built with MPI and OpenMP")
    endif()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
        # No fused multiply-add on the host or the GPU (as -ffp-contract=off).
        set(HEAT2D_OFFLOAD_DEFAULT "-mp=gpu;-gpu=cc86,nofma;-Minfo=mp")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(HEAT2D_OFFLOAD_DEFAULT "-fopenmp;--offload-arch=sm_86")
    else()
        set(HEAT2D_OFFLOAD_DEFAULT "")
        message(WARNING "HEAT2D_OFFLOAD: no GPU flags known for ${CMAKE_CXX_COMPILER_ID}; the target regions will run on the host")
    endif()
    set(HEAT2D_OFFLOAD_FLAGS "${HEAT2D_OFFLOAD_DEFAULT}" CACHE STRING "Offload compile and link flags for heat_parallel")
    target_compile_definitions(heat_parallel PRIVATE HEAT2D_USE_OFFLOAD=1)
    target_compile_options(heat_parallel PRIVATE ${HEAT2D_OFFLOAD_FLAGS})
    target_link_options(heat_parallel PRIVATE ${HEAT2D_OFFLOAD_FLAGS})
endif()

if(HEAT2D_BUILD_TESTS)
    enable_testing()
    foreach(HEAT2D_TEST velocity config grid io poisson flow affinity)
        add_executable(test_${HEAT2D_TEST} tests/test_${HEAT2D_TEST}.cpp)
        target_link_libraries(test_${HEAT2D_TEST} PRIVATE heat2d_common)
        add_test(NAME ${HEAT2D_TEST} COMMAND test_${HEAT2D_TEST})
        list(APPEND HEAT2D_TARGETS test_${HEAT2D_TEST})
    endforeach()

    add_executable(test_serial_solver tests/test_serial_solver.cpp src/serial/solver.cpp)
    target_link_libraries(test_serial_solver PRIVATE heat2d_common)
    add_test(NAME serial_solver COMMAND test_serial_solver)
    list(APPEND HEAT2D_TARGETS test_serial_solver)

    # The 3D pigment model against the 2D solver and its analytic references.
    add_executable(test_3d tests/test_3d.cpp src/serial/solver.cpp src/serial3d/solver.cpp)
    target_link_libraries(test_3d PRIVATE heat3d_common)
    add_test(NAME serial_3d COMMAND test_3d)
    list(APPEND HEAT2D_TARGETS test_3d)

    # End-to-end checks through the command-line interface (scripts/check_correctness.py).
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME correctness_serial
                 COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_correctness.py"
                         --serial-only --build-dir "${CMAKE_BINARY_DIR}"
                         --output-dir "${CMAKE_BINARY_DIR}/correctness"
                 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        add_test(NAME correctness_3d_serial
                 COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_correctness.py"
                         --experiment "${CMAKE_CURRENT_SOURCE_DIR}/experiments/correctness_3d.json"
                         --serial-only --build-dir "${CMAKE_BINARY_DIR}"
                         --output-dir "${CMAKE_BINARY_DIR}/correctness_3d"
                 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        # The 3D Poisson benchmark: second-order convergence and the exact discrete solution.
        add_test(NAME poisson3d_verify_serial
                 COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/poisson3d/run_experiments.py"
                         verify --serial-only --build-dir "${CMAKE_BINARY_DIR}"
                         --output "${CMAKE_BINARY_DIR}/poisson3d_verify"
                 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/poisson3d")
        # Serial vs parallel comparisons; needs up to 6 MPI ranks. Extra launcher
        # flags (e.g. --oversubscribe on a small machine) go in HEAT2D_MPIEXEC_FLAGS.
        if(TARGET heat_parallel AND MPIEXEC_EXECUTABLE)
            add_executable(test_halo tests/test_halo.cpp src/parallel/decomposition.cpp src/parallel/halo.cpp)
            target_link_libraries(test_halo PRIVATE heat2d_common MPI::MPI_CXX)
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                target_compile_options(test_halo PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
            endif()
            foreach(HEAT2D_RANKS 1 4 6)
                add_test(NAME halo_${HEAT2D_RANKS}_ranks
                         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${HEAT2D_RANKS} ${HEAT2D_MPIEXEC_FLAGS}
                                 $<TARGET_FILE:test_halo>)
            endforeach()

            # The 3D pigment halos (slabs and pencils, blocking and overlapping).
            add_executable(test_halo3d tests/test_halo3d.cpp src/parallel/decomposition.cpp src/parallel3d/halo3d.cpp)
            target_link_libraries(test_halo3d PRIVATE heat3d_common MPI::MPI_CXX)
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                target_compile_options(test_halo3d PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
            endif()
            foreach(HEAT2D_RANKS 1 4 6)
                add_test(NAME halo3d_${HEAT2D_RANKS}_ranks
                         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${HEAT2D_RANKS} ${HEAT2D_MPIEXEC_FLAGS}
                                 $<TARGET_FILE:test_halo3d>)
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
            add_test(NAME correctness_3d_parallel
                     COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_correctness.py"
                             --experiment "${CMAKE_CURRENT_SOURCE_DIR}/experiments/correctness_3d.json"
                             --build-dir "${CMAKE_BINARY_DIR}"
                             --output-dir "${CMAKE_BINARY_DIR}/correctness_3d"
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
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    # No fused multiply-add, as -ffp-contract=off above.
    foreach(HEAT2D_TARGET IN LISTS HEAT2D_TARGETS)
        target_compile_options(${HEAT2D_TARGET} PRIVATE -Mnofma)
    endforeach()
endif()
