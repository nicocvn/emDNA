# -------------------------------------------------------------------------- #
# tclap dependency CMake file
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# -------------------------------------------------------------------------- #


if(NOT TARGET tclap)
    add_library(tclap INTERFACE)
    target_include_directories(tclap INTERFACE
                               ${CMAKE_CURRENT_LIST_DIR}/tclap)
endif()
