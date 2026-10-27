cmake_minimum_required(VERSION 3.16)
project(ParallelHeat2D VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(HEAT2D_BUILD_SOLVERS "Build solvers after their source files have been implemented" OFF)
option(HEAT2D_ENABLE_MPI "Build the parallel MPI target" ON)
option(HEAT2D_ENABLE_OPENMP "Enable OpenMP in the parallel target" ON)

if(NOT HEAT2D_BUILD_SOLVERS)
    message(STATUS "Heat2D scaffold mode: no solver targets are generated.")
    message(STATUS "Implement the planned C++ files, then set HEAT2D_BUILD_SOLVERS=ON.")
    return()
endif()

set(HEAT2D_COMMON_SOURCES
    src/common/config.cpp
    src/common/velocity.cpp
    src/common/physics.cpp
    src/common/io.cpp
)
set(HEAT2D_SERIAL_SOURCES
    src/serial/main.cpp
    src/serial/solver.cpp
)
set(HEAT2D_PARALLEL_SOURCES
    src/parallel/main.cpp
    src/parallel/solver.cpp
    src/parallel/decomposition.cpp
    src/parallel/halo.cpp
)

set(HEAT2D_REQUIRED_SOURCES ${HEAT2D_COMMON_SOURCES} ${HEAT2D_SERIAL_SOURCES})
if(HEAT2D_ENABLE_MPI)
    list(APPEND HEAT2D_REQUIRED_SOURCES ${HEAT2D_PARALLEL_SOURCES})
endif()

set(HEAT2D_MISSING_SOURCES "")
foreach(HEAT2D_SOURCE IN LISTS HEAT2D_REQUIRED_SOURCES)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${HEAT2D_SOURCE}")
        list(APPEND HEAT2D_MISSING_SOURCES "${HEAT2D_SOURCE}")
    endif()
endforeach()
if(HEAT2D_MISSING_SOURCES)
    string(REPLACE ";" "\n  " HEAT2D_MISSING_TEXT "${HEAT2D_MISSING_SOURCES}")
    message(FATAL_ERROR
        "Solver sources have not been implemented yet:\n  ${HEAT2D_MISSING_TEXT}\n"
        "Use HEAT2D_BUILD_SOLVERS=OFF to configure the scaffold."
    )
endif()

add_library(heat2d_common STATIC ${HEAT2D_COMMON_SOURCES})
target_include_directories(heat2d_common PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_executable(heat_serial ${HEAT2D_SERIAL_SOURCES})
target_link_libraries(heat_serial PRIVATE heat2d_common)

set(HEAT2D_TARGETS heat2d_common heat_serial)
if(HEAT2D_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    add_executable(heat_parallel ${HEAT2D_PARALLEL_SOURCES})
    target_link_libraries(heat_parallel PRIVATE heat2d_common MPI::MPI_CXX)

    if(HEAT2D_ENABLE_OPENMP)
        find_package(OpenMP REQUIRED COMPONENTS CXX)
        target_link_libraries(heat_parallel PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(heat_parallel PRIVATE HEAT2D_USE_OPENMP=1)
    endif()
    list(APPEND HEAT2D_TARGETS heat_parallel)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(HEAT2D_TARGET IN LISTS HEAT2D_TARGETS)
        target_compile_options(${HEAT2D_TARGET} PRIVATE -Wall -Wextra -Wpedantic)
    endforeach()
endif()
