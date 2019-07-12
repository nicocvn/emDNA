# DNASimTests CMake


# --- CMake gobal directives ---

cmake_minimum_required (VERSION 2.6)


# --- DNASim project ---

project (DNASimTests CXX)
include (${CMAKE_SOURCE_DIR}/DNASimCommon.cmake)


# --- build type options ---
# disable googletest threading
add_definitions(-DGTEST_HAS_PTHREAD=0 -DGTEST_USE_OWN_TR1_TUPLE=1)


# --- source tree directives --- #

# add all source directories
set (SourceDirectories
	.
	)
set (GtestDirectory
	gtest)

# cpp files list
append_files (DNASimTestsCppSources "*cpp" ${SourceDirectories})
append_files (GtestCppSources "*cc" ${GtestDirectory})

# include directories
include_directories (BEFORE .)
include_directories (BEFORE ../src)
include_directories (BEFORE ../src/dna)
include_directories (BEFORE ../src/file_io)
include_directories (BEFORE ../src/geometry)
include_directories (BEFORE ../src/maths)
include_directories (BEFORE ../src/prn-generators)
include_directories (BEFORE ../src/serialization)
include_directories (BEFORE ../src/simulations)
include_directories (BEFORE ../deps/)

# ODE support option
option (DNASIM_ODE_COLLISION "Include ODE collision detection." OFF)
if (DNASIM_ODE_COLLISION)

	# macro and headers
	add_definitions(-DWITH_ODE_COLLISION)
	include_directories (BEFORE deps/opende/include)
	
endif (DNASIM_ODE_COLLISION)


# --- target directives ---

link_directories (${CMAKE_BINARY_DIR})
add_executable (DNASimTests ${DNASimTestsCppSources} ${GtestCppSources})
target_link_libraries (DNASimTests DNASim)

