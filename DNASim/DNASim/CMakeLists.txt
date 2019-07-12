# DNASim CMake


# --- CMake gobal directives ---

cmake_minimum_required (VERSION 2.6)


# --- DNASim project ---

project (DNASim C CXX)
include (DNASimCommon.cmake)
set (DNASIM_MAJOR_VERSION "0")
set (DNASIM_MINOR_VERSION "0")
set (DNASIM_BUGFIX_VERSION "0")
set (DNASIM_VERSION_NUMBER
    "${DNASIM_MAJOR_VERSION}.${DNASIM_MINOR_VERSION}.${DNASIM_BUGFIX_VERSION}"
    )
set (DNASIM_VERSION "${DNASIM_VERSION_NUMBER}")


# --- build type options ---

if (NOT CMAKE_BUILD_TYPE)

	# compiler
	# force to discard any user input
	set (COMPILER ${CMAKE_CXX_COMPILER})
	mark_as_advanced (CLEAR CMAKE_CXX_COMPILER)
	set (CMAKE_CXX_COMPILER ${COMPILER}
		CACHE
		STRING
		"C++ compiler"
		FORCE)

	# set build type as MacOS
	# force to discard any user input
	set (CMAKE_BUILD_TYPE MACOS
		CACHE
		STRING
    	"MacOS build settings."
    	FORCE)

	# set MacOS c++ compiler flags
    set (MACOS_FLAGS "-O2 -march=native -msse")
    set (MACOS_FLAGS "${MACOS_FLAGS} -std=gnu++11 -stdlib=libc++")
	set (MACOS_FLAGS "${MACOS_FLAGS} -m64 -funroll-loops")
    set (MACOS_FLAGS "${MACOS_FLAGS} -msse3 -mssse3 -msse4.1 -msse4.2")
	set (MACOS_FLAGS "${MACOS_FLAGS}  -ffast-math -fstrict-aliasing")
	set (MACOS_FLAGS "${MACOS_FLAGS} -Wno-deprecated -Weffc++")
	set (CMAKE_CXX_FLAGS_MACOS ${MACOS_FLAGS}
		CACHE
		FILEPATH
		"C++ MacOS compiler flags.")

	# set Unix cluster c++ compiler flags
	set (UNIX_FLAGS "-O2 -march=native -msse -mfpmath=sse")
	set (UNIX_FLAGS "${UNIX_FLAGS} -std=gnu++11")
	set (UNIX_FLAGS "${UNIX_FLAGS} -m64 -funroll-loops")
    set (UNIX_FLAGS "${UNIX_FLAGS} -msse3 -mssse3 -ffast-math")
	set (UNIX_FLAGS "${UNIX_FLAGS} -fstrict-aliasing")
	set (UNIX_FLAGS "${UNIX_FLAGS} -Wno-deprecated")
	set (CMAKE_CXX_FLAGS_UNIX ${UNIX_FLAGS}
		CACHE
		FILEPATH
		"C++ Unix compiler flags (optimized).")

	# OS X architectures
	if (APPLE)
		set (CMAKE_OSX_ARCHITECTURES "x86_64"
			CACHE
			STRING
			"Architecture to build against (OS X only)."
			FORCE)
        set (CMAKE_MACOSX_RPATH 1)
	endif (APPLE)

endif (NOT CMAKE_BUILD_TYPE)

if (APPLE)
    set (CMAKE_MACOSX_RPATH 1)
endif (APPLE)


# --- source tree directives --- #

# add all source directories
set (SourceDirectories
	src/dna
	src/file_io
	src/geometry
	src/maths
    src/pdb
	src/prn-generators
	src/serialization
	src/simulations
	)

# headers and cpp files list
append_files (DNASimHeaders "h" ${SourceDirectories})
append_files (DNASimCppSources "cpp" ${SourceDirectories})

# include directories
include_directories (BEFORE src)
include_directories (BEFORE src/dna)
include_directories (BEFORE src/file_io)
include_directories (BEFORE src/geometry)
include_directories (BEFORE src/maths)
include_directories (BEFORE src/pdb)
include_directories (BEFORE src/prn-generators)
include_directories (BEFORE src/serialization)
include_directories (BEFORE src/simulations)
include_directories (BEFORE src/dna_force_field_packager)
include_directories (BEFORE deps/)

# ODE support option
option (DNASIM_ODE_COLLISION "Include ODE collision detection." OFF)
if (DNASIM_ODE_COLLISION)

    set(ODE_TAR "${CMAKE_SOURCE_DIR}/deps/ode-0.13/ode-0.13.tar.gz")
    set(ODE_SRC_DIR "ode-0.13")

	# ODE static library target
	add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/ode_shared/lib/libode.a
        COMMAND rm -Rf ${ODE_SRC_DIR}
        COMMAND tar xf ${ODE_TAR}
		COMMAND cp ${CMAKE_SOURCE_DIR}/deps/opende_compile_script.sh .
		COMMAND sh opende_compile_script.sh shared
        COMMAND rm -Rf ${ODE_SRC_DIR}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Building ODE static library (fPIC)."
		)
    add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/ode_static/lib/libode.a
        COMMAND rm -Rf ${ODE_SRC_DIR}
        COMMAND tar xf ${ODE_TAR}
        COMMAND cp ${CMAKE_SOURCE_DIR}/deps/opende_compile_script.sh .
        COMMAND sh opende_compile_script.sh static
        COMMAND rm -Rf ${ODE_SRC_DIR}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Building ODE static library."
    )

    # macro and headers
    add_definitions(-DWITH_ODE_COLLISION)
    include_directories (BEFORE ${CMAKE_BINARY_DIR}/ode_shared/include)
    include_directories (BEFORE ${CMAKE_BINARY_DIR}/ode_static/include)

endif (DNASIM_ODE_COLLISION)


# --- target directives ---

# ODE support
if (DNASIM_ODE_COLLISION)
    add_library (DNASim SHARED
        ${DNASimHeaders}
        ${DNASimCppSources}
        ${CMAKE_BINARY_DIR}/ode_shared/lib/libode.a
        )
    add_library (DNASimStatic STATIC
        ${DNASimHeaders}
        ${DNASimCppSources}
        ${CMAKE_BINARY_DIR}/ode_static/lib/libode.a
        )
    target_link_libraries (DNASim
        ${CMAKE_BINARY_DIR}/ode_shared/lib/libode.a
        )
    target_link_libraries (DNASimStatic
        ${CMAKE_BINARY_DIR}/ode_shared/lib/libode.a
        )

    # static objects
    add_custom_command(TARGET DNASimStatic
        POST_BUILD
        COMMAND ar x ode_static/lib/libode.a
        COMMAND ar rs libDNASim_static.a *.o
        COMMAND rm -Rf *.o __.SYMDEF
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Assembling object files."
        )

# no ODE
else (DNASIM_ODE_COLLISION)
    add_library (DNASim SHARED
        ${DNASimHeaders}
        ${DNASimCppSources}
        )
    add_library (DNASimStatic STATIC
        ${DNASimHeaders}
        ${DNASimCppSources}
        )
endif (DNASIM_ODE_COLLISION)

# add _static to the static library name
# set visibility as hidden
set_target_properties(DNASimStatic
    PROPERTIES
    OUTPUT_NAME DNASim_static
    COMPILER_FLAGS "-fno-visibility-inlines-hidden -fvisibility=default"
    )

# add fPIC flag for shared lib objects
# set visibility as hidden
# rpath for shared lib
set_target_properties(DNASim
    PROPERTIES
    COMPILER_FLAGS "-fPIC -fno-visibility-inlines-hidden -fvisibility=default"
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )

# ffield packager
add_executable(ffield_packager
    src/dna_force_field_packager/ffield_packager.cpp)
target_link_libraries(ffield_packager DNASimStatic)


# --- testing ---

option (DNASIM_TEST "Builds DNASim unit tests." OFF)
if (DNASIM_TEST)
	add_subdirectory(tests)
endif (DNASIM_TEST)


# --- install directives ---

# library install
set (DNASimLibPermissions OWNER_READ;OWNER_WRITE;OWNER_EXECUTE)
set (DNASimLibPermissions ${DNASimLibPermissions} GROUP_READ;GROUP_EXECUTE)
set (DNASimLibPermissions ${DNASimLibPermissions} WORLD_READ;WORLD_EXECUTE)
install (TARGETS DNASim
	DESTINATION lib
	PERMISSIONS ${DNASimLibPermissions}
	)
install (TARGETS DNASimStatic
    DESTINATION lib
    PERMISSIONS ${DNASimLibPermissions}
    )

# headers install
set (HeadersDirectories
	src
	src/dna
	src/file_io
	src/geometry
	src/maths
    src/pdb
	src/prn-generators
    src/serialization
	src/simulations
	)
append_files (DNASimLibHeaders "h" ${HeadersDirectories})
set (DNASimHeadersPermissions OWNER_READ;OWNER_WRITE)
set (DNASimHeadersPermissions ${DNASimHeadersPermissions} GROUP_READ)
set (DNASimHeadersPermissions ${DNASimHeadersPermissions} WORLD_READ)
install (FILES ${DNASimLibHeaders}
	DESTINATION include/DNASim
    PERMISSIONS ${DNASimHeadersPermissions}
	)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/deps/Eigen
	DESTINATION include/DNASim
    )
install (DIRECTORY ${CMAKE_SOURCE_DIR}/deps/tclap
    DESTINATION include/DNASim
    )
install (DIRECTORY ${CMAKE_SOURCE_DIR}/deps/cereal
    DESTINATION include/DNASim
    )
install (DIRECTORY ${CMAKE_SOURCE_DIR}/deps/nanoflann
    DESTINATION include/DNASim
    )
if (DNASIM_ODE_COLLISION)
	install (DIRECTORY ${CMAKE_BINARY_DIR}/ode_shared/include/ode
		DESTINATION include/DNASim
        )
endif (DNASIM_ODE_COLLISION)

# dna force field packager
set (toolPermissions OWNER_READ;OWNER_WRITE;OWNER_EXECUTE)
set (toolPermissions ${toolPermissions} GROUP_READ;GROUP_EXECUTE)
set (toolPermissions ${toolPermissions} WORLD_READ;WORLD_EXECUTE)
install (TARGETS ffield_packager
    DESTINATION bin
    PERMISSIONS ${toolPermissions}
    )


# --- status messages ---
message (STATUS "")
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message (STATUS
		"C++ compiler (CMAKE_CXX_COMPILER): ${CMAKE_CXX_COMPILER}")
message (STATUS
		"C++ compiler flags (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}): "
		"${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
message (STATUS "Installation path (CMAKE_INSTALL_PREFIX): "
		"${CMAKE_INSTALL_PREFIX}")
if (APPLE)
	message (STATUS "OS X architecture to build (CMAKE_OSX_ARCHITECTURES): "
			"${CMAKE_OSX_ARCHITECTURES}")
	message (STATUS "OS X system target (CMAKE_OSX_DEPLOYMENT_TARGET): "
			"${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif (APPLE)
message (STATUS "Include ODE collision detection (DNASIM_ODE_COLLISION): "
	"${DNASIM_ODE_COLLISION}")
message (STATUS "Builds DNASim unit tests (DNASIM_TEST): "
	"${DNASIM_TEST}")
message (STATUS "")
if (DNASIM_TEST)
	message (STATUS "Run ${CMAKE_BINARY_DIR}/tests/DNASimTests to test the "
		"library\n")
endif (DNASIM_TEST)
