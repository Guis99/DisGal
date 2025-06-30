cmake_minimum_required(VERSION 3.15)
project(oneDsrc)

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(../include)
include_directories(../../include)
include_directories(../../Dependencies)
include_directories(../../Utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

set(CURRENT_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build)

cmake_print_variables(CMAKE_CURRENT_LIST_DIR)
cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

set(DEPENDENCIES_ROOT ${CURRENT_SRC_DIR}/../../Dependencies)

cmake_print_variables(DEPENDENCIES_ROOT)

# Define source directories
set(MATHPARSER_SRC ${DEPENDENCIES_ROOT}/MathParser/src)
set(QTM_SRC ${DEPENDENCIES_ROOT}/QTM/src)
set(UTILS_SRC ${CURRENT_SRC_DIR}/../../Utils)

file(GLOB MATHPARSER_SRC_FILES "${MATHPARSER_SRC}/*.cpp")
file(GLOB QTM_SRC_FILES "${QTM_SRC}/*.cpp")
file(GLOB UTILS_SRC_FILES "${UTILS_SRC}/*.cpp")

cmake_print_variables(CURRENT_SRC_DIR)

add_library(oneDsrc STATIC 
    ${CURRENT_SRC_DIR}/common.cpp)

set_target_properties(oneDsrc PROPERTIES OBJECT_OUTPUTS "${BUILD_DIR}/")

add_library(mathparser OBJECT ${MATHPARSER_SRC_FILES})
set_target_properties(mathparser PROPERTIES OBJECT_OUTPUTS "${BUILD_DIR}/")

add_library(qtm OBJECT ${QTM_SRC_FILES})
set_target_properties(qtm PROPERTIES OBJECT_OUTPUTS "${BUILD_DIR}/")

add_library(utils OBJECT ${UTILS_SRC_FILES})
set_target_properties(utils PROPERTIES OBJECT_OUTPUTS "${BUILD_DIR}/")

target_link_libraries(oneDsrc 
    PRIVATE 
        mathparser 
        qtm 
        utils)