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

endforeach()

# Static library.
add_library(DNASim STATIC ${DNASim_SOURCES})
target_link_libraries(DNASim PUBLIC cereal eigen nanoflann)
target_include_directories(DNASim PUBLIC ${DNASim_HEADERS_DIRS})

