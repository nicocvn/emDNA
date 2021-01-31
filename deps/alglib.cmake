# -------------------------------------------------------------------------- #
# alglib dependency CMake file
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# -------------------------------------------------------------------------- #


if(NOT TARGET alglib)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/alglib/)
endif()
