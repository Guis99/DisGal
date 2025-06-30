cmake_minimum_required(VERSION 3.15)
project(PoissonSRC)

include(CMakePrintHelpers)
include(../../src/CMakeLists.txt)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(../../include)
include_directories(../../../include)
include_directories(../../../Dependencies)
include_directories(../../../Utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")

set(MAIN_SRC_FILE main.cpp)
add_executable(main ${MAIN_SRC_FILE})
target_link_libraries(main PRIVATE oneDsrc)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})