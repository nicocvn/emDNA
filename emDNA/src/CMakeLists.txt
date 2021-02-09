# ---------------------------------------------------------------------------- #
# emDNA library
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# --- List of project subdirectories ---

# List of subdirectories in project tree.
set(source_tree
    DNA
    lego
    minim
    utils)


# --- Static library target ---

# Add all subdirectories and grab all sources and headers directories.
set(emDNA_SOURCES emDNA.h)
set(emDNA_HEADERS_DIRS ${CMAKE_CURRENT_LIST_DIR})
foreach(dir IN LISTS source_tree)

    # Load subdirectory.
    add_subdirectory(${dir})

    # Append source files and headers directories.
    list(APPEND emDNA_SOURCES ${${dir}_SOURCES})

endforeach()

# Static library.
add_library(emDNA STATIC ${emDNA_SOURCES})
target_link_libraries(emDNA PUBLIC DNASim alglib)
target_include_directories(emDNA PUBLIC ${emDNA_HEADERS_DIRS})
