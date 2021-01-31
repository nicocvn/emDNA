# ---------------------------------------------------------------------------- #
# DNASim library
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# --- List of project subdirectories ---

# List of subdirectories in project tree.
set(source_tree
    dna
    file_io
    geometry
    maths
    pdb
    prn-generators
    serialization
    simulations)


# --- Static library target ---

# Add all subdirectories and grab all sources and headers directories.
set(DNASim_SOURCES)
set(DNASim_HEADERS_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(dir IN LISTS source_tree)

    # Load subdirectory.
    add_subdirectory(${dir})

    # Append source files and headers directories.
    list(APPEND DNASim_SOURCES ${${dir}_SOURCES})
#    list(APPEND sdk_headers_dirs ${${dir}_headers_dirs})

endforeach()

# Static library.
add_library(DNASim STATIC ${DNASim_SOURCES})
target_link_libraries(DNASim PUBLIC cereal eigen nanoflann)
# todo: edit once Eigen is better managed as a dependency
target_include_directories(DNASim PUBLIC ${DNASim_HEADERS_DIRS} ../deps)


# --- Install directives ---

# Destination variables.
set(include_dest "include/DNASim")
set(main_lib_dest "lib/DNASim")
set(resource_dest "share/DNASim")

# Include directories interfaces for libraries.
set_property(TARGET DNASim
             PROPERTY
             INTERFACE_INCLUDE_DIRECTORIES
             $<BUILD_INTERFACE:${DNASim_HEADERS_DIRS}>
             $<INSTALL_INTERFACE:${include_dest}>)

# Headers.
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ DESTINATION ${include_dest}
        FILES_MATCHING PATTERN "*.h")

# Library.
install(TARGETS DNASim
        CONFIGURATIONS ${CMAKE_BUILD_TYPE}
        DESTINATION ${main_lib_dest}/${CMAKE_BUILD_TYPE}
        EXPORT emDNA)

## CMake export file.
#install(EXPORT emDNA
#        NAMESPACE emDNA::
#        DESTINATION ${main_lib_dest}/${CMAKE_BUILD_TYPE})
