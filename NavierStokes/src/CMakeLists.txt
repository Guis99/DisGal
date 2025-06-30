cmake_minimum_required(VERSION 3.15)
project(PoissonSRC)

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

set(CURRENT_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build)

set(DEPENDENCIES_ROOT ${CURRENT_SRC_DIR}/../../Dependencies)
set(UTILS_ROOT ${CURRENT_SRC_DIR}/../../Utils)

# Define source directories
set(MATHPARSER_SRC ${DEPENDENCIES_ROOT}/MathParser/src)
set(QTM_SRC ${DEPENDENCIES_ROOT}/QTM/src)
set(UTILS_SRC ${CURRENT_SRC_DIR}/../../Utils)

file(GLOB MATHPARSER_SRC_FILES "${MATHPARSER_SRC}/*.cpp")
file(GLOB QTM_SRC_FILES "${QTM_SRC}/*.cpp")
file(GLOB UTILS_SRC_FILES "${UTILS_SRC}/*.cpp")

add_library(poissonsrc STATIC 
    ${CURRENT_SRC_DIR}/NSMatrixAssemblydG.cpp)

target_include_directories(poissonsrc 
    PRIVATE
        ${CURRENT_SRC_DIR}/../include
        ${UTILS_ROOT}/include
        ${DEPENDENCIES_ROOT})

target_sources(poissonsrc 
    PRIVATE
        ${MATHPARSER_SRC_FILES}
        ${QTM_SRC_FILES}
        ${UTILS_SRC_FILES})

set_target_properties(poissonsrc PROPERTIES OBJECT_OUTPUTS "${BUILD_DIR}/")