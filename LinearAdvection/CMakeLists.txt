cmake_minimum_required(VERSION 3.15)
project(LinearAdvection)

# project root = DisGal/LinearAdvection

include(CMakePrintHelpers)

# set(CMAKE_C_COMPILER /usr/bin/clang)
# set(CMAKE_CXX_COMPILER /usr/bin/clang++)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Append compiler flags instead of overwriting
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DVERBOSE -std=c++17 -stdlib=libc++ -nostdinc++ -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk/usr/include/c++/v1")

# std=c++17 -stdlib=libc++ -I/Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk/usr/include/c++/v1

# Ensure standard library paths are included
if(APPLE)
    # include_directories(/Library/Developer/CommandLineTools/usr/include/c++/v1)
elseif(UNIX)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/c++/$(gcc -dumpversion)")
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
set(DEPENDENCIES_ROOT ${PROJECT_ROOT}/../Dependencies)

cmake_print_variables(DEPENDENCIES_ROOT)

# Define build directories
set(LINEAR_ADVECTION_BUILD_DIR ${PROJECT_ROOT}/build)
set(MATHPARSER_BUILD_DIR ${DEPENDENCIES_ROOT}/MathParser/build)
set(UTILS_BUILD_DIR ${PROJECT_ROOT}/../Utils/build)

# Ensure each build directory exists
file(MAKE_DIRECTORY ${LINEAR_ADVECTION_BUILD_DIR})
file(MAKE_DIRECTORY ${MATHPARSER_BUILD_DIR})
file(MAKE_DIRECTORY ${UTILS_BUILD_DIR})

# Source directories
set(LINEAR_ADVECTION_SRC ${PROJECT_ROOT})
set(MATHPARSER_SRC ${DEPENDENCIES_ROOT}/MathParser/src)
set(UTILS_SRC ${LINEAR_ADVECTION_SRC}/../Utils)

# Include directories
include_directories(${LINEAR_ADVECTION_SRC})
include_directories(${DEPENDENCIES_ROOT}/MathParser/include)
include_directories(${UTILS_SRC})

# Source files
set(COMMON_SRC ${LINEAR_ADVECTION_SRC}/common.cpp)
set(LINEAR_ADVECTION_MAIN_SRC ${LINEAR_ADVECTION_SRC}/1D.cpp)  # Change to 1D.cpp to compile that instead
file(GLOB MATHPARSER_SRC_FILES "${MATHPARSER_SRC}/*.cpp")
file(GLOB UTILS_SRC_FILES "${UTILS_SRC}/*.cpp")

cmake_print_variables(UTILS_SRC)
cmake_print_variables(UTILS_SRC_FILES)

# Compile each source file separately into object files in respective build directories
add_library(common OBJECT ${COMMON_SRC})
set_target_properties(common PROPERTIES OBJECT_OUTPUTS "${LINEAR_ADVECTION_BUILD_DIR}/common.o")

add_library(mathparser OBJECT ${MATHPARSER_SRC_FILES})
set_target_properties(mathparser PROPERTIES OBJECT_OUTPUTS "${MATHPARSER_BUILD_DIR}/")

add_library(utils OBJECT ${UTILS_SRC_FILES})
set_target_properties(utils PROPERTIES OBJECT_OUTPUTS "${UTILS_BUILD_DIR}/")

# Compile main source file into an object file
set(MAIN_OBJECT "${LINEAR_ADVECTION_BUILD_DIR}/Adv1D.o")

cmake_print_variables(LINEAR_ADVECTION_BUILD_DIR)
cmake_print_variables(MAIN_OBJECT)
add_custom_command(
    OUTPUT ${MAIN_OBJECT}
    COMMAND ${CMAKE_CXX_COMPILER} -std=c++17 -c ${LINEAR_ADVECTION_MAIN_SRC} -o ${MAIN_OBJECT} -O3 -DVERBOSE
    DEPENDS ${LINEAR_ADVECTION_MAIN_SRC}
)

# Final executable
add_executable(Adv1D ${MAIN_OBJECT} $<TARGET_OBJECTS:common> $<TARGET_OBJECTS:mathparser> $<TARGET_OBJECTS:utils>)

# Set output location
set_target_properties(Adv1D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${LINEAR_ADVECTION_BUILD_DIR})

# Link dependencies
target_link_libraries(Adv1D PRIVATE common mathparser utils)

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
