# ---------------------------------------------------------------------------- #
# emDNA CMake project
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# --- Project setup ---

# CMake minimal version.
cmake_minimum_required(VERSION 3.13)

# Support for ccache.
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM}
        CACHE STRING "(emDNA) C compiler ccache program" FORCE)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM}
        CACHE STRING "(emDNA) C++ compiler ccache program" FORCE)
    message("C/C++ compiler ccache: ${CCACHE_PROGRAM}")
endif()

# Define project.
project(emDNA C CXX)

# Enforce C++14 standard.
set(CMAKE_CXX_STANDARD 14)

# Enforce C99 standard.
set(CMAKE_C_STANDARD 99)


# --- Project components ---

# Dependencies.
include(deps/alglib.cmake)
include(deps/cereal.cmake)
include(deps/eigen.cmake)
include(deps/googletest.cmake)
include(deps/nanoflann.cmake)
include(deps/tclap.cmake)

# DNASim
add_subdirectory(DNASim/)

# emDNA
add_subdirectory(emDNA/)