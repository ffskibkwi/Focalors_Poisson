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

add_subdirectory(solver)
add_subdirectory(case)
