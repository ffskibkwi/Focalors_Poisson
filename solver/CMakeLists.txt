file(GLOB_RECURSE SOLVER_HEADER_FILES "*.h" "*.hpp")
file(GLOB_RECURSE SOLVER_SOURCE_FILES "*.cpp")

add_library(Solver ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES})
target_include_directories(Solver PUBLIC ${CMAKE_SOURCE_DIR}/solver)
target_precompile_headers(Solver PUBLIC pch.h)

# Avoid following error:
# Policy CMP0074 is not set: find_package uses <PackageName>_ROOT variables.
cmake_policy(SET CMP0074 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(FFTW REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI) # MPI_C may not found, so don't use REQUIRED

# parallel

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