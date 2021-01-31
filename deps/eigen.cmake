# -------------------------------------------------------------------------- #
# Eigen dependency CMake file
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# -------------------------------------------------------------------------- #


# Remarks:
# we only use a very limited scope of Eigen code base so rather than going
# full blast on the Eigen project we only create a CMake interface library
# which is sufficient for our needs.


if(NOT TARGET eigen)
    add_library(eigen INTERFACE)
    target_include_directories(eigen INTERFACE
                               ${CMAKE_CURRENT_LIST_DIR}/eigen)
endif()
