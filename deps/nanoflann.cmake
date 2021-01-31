# -------------------------------------------------------------------------- #
# nanoflann dependency CMake file
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# -------------------------------------------------------------------------- #


if(NOT TARGET nanoflann)
    add_library(nanoflann INTERFACE)
    target_include_directories(nanoflann INTERFACE
                               ${CMAKE_CURRENT_LIST_DIR}/nanoflann/include)
endif()
