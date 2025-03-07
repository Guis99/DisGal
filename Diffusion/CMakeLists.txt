cmake_minimum_required(VERSION 3.15)
project(Diffusion)

# project root = DisGal/Diffusion

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    set(CMAKE_C_COMPILER /usr/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-value -Wno-switch -O3 -DVERBOSE -std=c++17 -stdlib=libc++")
elseif(UNIX)
    set(CMAKE_C_COMPILER /usr/bin/gcc)
    set(CMAKE_CXX_COMPILER /usr/bin/g++)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wno-unused-value -Wno-switch -O3 -DVERBOSE -DMULTITHREAD -std=c++17 -stdlib=libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/c++/$(gcc -dumpversion)")
endif()

cmake_print_variables(CMAKE_CXX_FLAGS)

# Check if <vector> is available
include(CheckIncludeFileCXX)
check_include_file_cxx(vector HAVE_VECTOR)
if(NOT HAVE_VECTOR)
    message(FATAL_ERROR "The standard C++ <vector> header is missing! Check your compiler setup.")
endif()


# Set root paths
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
cmake_print_variables(PROJECT_ROOT)
set(DEPENDENCIES_ROOT ${PROJECT_ROOT}/../Dependencies)

cmake_print_variables(DEPENDENCIES_ROOT)

# Define build directories
set(DG_POISSON_BUILD_DIR ${PROJECT_ROOT}/build)
set(MATHPARSER_BUILD_DIR ${DEPENDENCIES_ROOT}/MathParser/build)
set(QTM_BUILD_DIR ${DEPENDENCIES_ROOT}/QTM/build)
set(UTILS_BUILD_DIR ${PROJECT_ROOT}/../Utils/build)

# Ensure each build directory exists
file(MAKE_DIRECTORY ${DG_POISSON_BUILD_DIR})
file(MAKE_DIRECTORY ${MATHPARSER_BUILD_DIR})
file(MAKE_DIRECTORY ${QTM_BUILD_DIR})
file(MAKE_DIRECTORY ${UTILS_BUILD_DIR})

# Source directories
set(DG_POISSON_SRC ${PROJECT_ROOT}/src)
set(MATHPARSER_SRC ${DEPENDENCIES_ROOT}/MathParser/src)
set(QTM_SRC ${DEPENDENCIES_ROOT}/QTM/src)
set(UTILS_SRC ${PROJECT_ROOT}/../Utils)

# Include directories
include_directories(${PROJECT_ROOT}/include)
include_directories(${DEPENDENCIES_ROOT}/MathParser/include)
include_directories(${DEPENDENCIES_ROOT}/QTM/include)
include_directories(${UTILS_SRC})

# Source files
set(DG_POISSON_MAIN_SRC ${PROJECT_ROOT}/mainTestPoissondG.cpp)  # Change to mainTest.cpp to compile that instead
file(GLOB DG_POISSON_SRC_FILES "${DG_POISSON_SRC}/*.cpp")
file(GLOB MATHPARSER_SRC_FILES "${MATHPARSER_SRC}/*.cpp")
file(GLOB QTM_SRC_FILES "${QTM_SRC}/*.cpp")
file(GLOB UTILS_SRC_FILES "${UTILS_SRC}/*.cpp")

cmake_print_variables(UTILS_SRC)
cmake_print_variables(UTILS_SRC_FILES)

# Compile each source file separately into object files in respective build directories

add_library(diffusiondg OBJECT ${DG_POISSON_SRC_FILES})
set_target_properties(diffusiondg PROPERTIES OBJECT_OUTPUTS "${DG_POISSON_BUILD_DIR}/")

add_library(mathparser OBJECT ${MATHPARSER_SRC_FILES})
set_target_properties(mathparser PROPERTIES OBJECT_OUTPUTS "${MATHPARSER_BUILD_DIR}/")

add_library(qtm OBJECT ${QTM_SRC_FILES})
set_target_properties(qtm PROPERTIES OBJECT_OUTPUTS "${QTM_BUILD_DIR}/")

add_library(utils OBJECT ${UTILS_SRC_FILES})
set_target_properties(utils PROPERTIES OBJECT_OUTPUTS "${UTILS_BUILD_DIR}/")

# Compile main source file into an object file
set(MAIN_OBJECT "${DG_POISSON_BUILD_DIR}/diffDG.o")

cmake_print_variables(DG_POISSON_BUILD_DIR)
cmake_print_variables(MAIN_OBJECT)
add_custom_command(
    OUTPUT ${MAIN_OBJECT}
    COMMAND ${CMAKE_CXX_COMPILER} -std=c++17 -c ${DG_POISSON_MAIN_SRC} -o ${MAIN_OBJECT} -O3 -DVERBOSE
    DEPENDS ${DG_POISSON_MAIN_SRC}
)

# Final executable
add_executable(diffDG ${MAIN_OBJECT} $<TARGET_OBJECTS:diffusiondg> $<TARGET_OBJECTS:mathparser> $<TARGET_OBJECTS:qtm> $<TARGET_OBJECTS:utils>)

# Set output location
set_target_properties(diffDG PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${DG_POISSON_BUILD_DIR})

# Link dependencies
target_link_libraries(diffDG PRIVATE diffusiondg mathparser qtm utils)

# Platform-specific settings
if(APPLE)
    message(STATUS "Building for macOS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
elseif(UNIX)
    message(STATUS "Building for Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
elseif(WIN32)
    message(STATUS "Building for Windows")
endif()
