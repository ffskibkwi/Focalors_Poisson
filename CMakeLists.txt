cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CONFIGURATION_TYPES
    Debug Hybrid PureMPI
    CACHE STRING "Available build configurations" FORCE)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
if(CMAKE_BUILD_TYPE STREQUAL "Hybrid")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -DHYBRID -DNDEBUG")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "PureMPI")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -DPUREMPI -DNDEBUG")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ---------------------- Solver library (moved from solver/CMakeLists.txt) ----------------------
file(GLOB_RECURSE SOLVER_HEADER_FILES
    "${CMAKE_SOURCE_DIR}/base/*.h" "${CMAKE_SOURCE_DIR}/base/*.hpp"
    "${CMAKE_SOURCE_DIR}/io/*.h"   "${CMAKE_SOURCE_DIR}/io/*.hpp"
    "${CMAKE_SOURCE_DIR}/pe/*.h"   "${CMAKE_SOURCE_DIR}/pe/*.hpp"
    "${CMAKE_SOURCE_DIR}/pe_mpi/*.h" "${CMAKE_SOURCE_DIR}/pe_mpi/*.hpp"
)
file(GLOB_RECURSE SOLVER_SOURCE_FILES
    "${CMAKE_SOURCE_DIR}/base/*.cpp"
    "${CMAKE_SOURCE_DIR}/io/*.cpp"
    "${CMAKE_SOURCE_DIR}/pe/*.cpp"
    "${CMAKE_SOURCE_DIR}/pe_mpi/*.cpp"
)

# 排除 pe/poisson_old 和 pe/concat_old 目录下所有 cpp 文件
file(GLOB_RECURSE POISSON_OLD_SOURCE_FILES "${CMAKE_SOURCE_DIR}/pe/poisson_old/*.cpp")
file(GLOB_RECURSE CONCAT_OLD_SOURCE_FILES "${CMAKE_SOURCE_DIR}/pe/concat_old/*.cpp")
list(REMOVE_ITEM SOLVER_SOURCE_FILES ${POISSON_OLD_SOURCE_FILES})
list(REMOVE_ITEM SOLVER_SOURCE_FILES ${CONCAT_OLD_SOURCE_FILES})

add_library(Solver ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES})
target_include_directories(Solver PUBLIC ${CMAKE_SOURCE_DIR})
target_precompile_headers(Solver PUBLIC base/pch.h)

# Avoid following error:
# Policy CMP0074 is not set: find_package uses <PackageName>_ROOT variables.
cmake_policy(SET CMP0074 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(FFTW REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI) # MPI_C may not found, so don't use REQUIRED

target_include_directories(Solver PUBLIC
    ${FFTW_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_DIRS})
target_link_libraries(Solver PUBLIC
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,19.14>>:stdc++fs>
)
target_link_libraries(Solver PUBLIC
    ${FFTW_LIBRARIES}
    OpenMP::OpenMP_CXX
    MPI::MPI_CXX)

# ------------------------------------------------------------------------------------------------

add_subdirectory(case)
