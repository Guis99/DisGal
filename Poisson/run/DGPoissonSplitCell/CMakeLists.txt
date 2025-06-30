cmake_minimum_required(VERSION 3.15)
project(PoissonSRC_splitcell)

include(CMakePrintHelpers)
include(../../src/CMakeLists.txt)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR}/../..)

set(MAIN_SRC_FILE main.cpp)
add_executable(main ${MAIN_SRC_FILE})
target_link_libraries(main PRIVATE poissonsrc)
target_include_directories(main 
    PRIVATE
        ${PROJECT_ROOT}/include
        ${PROJECT_ROOT}/../Utils/include
        ${DEPENDENCIES_ROOT})