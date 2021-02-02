# -------------------------------------------------------------------------- #
# Cereal dependency CMake file
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# -------------------------------------------------------------------------- #


if(NOT TARGET cereal)

    # Options:
    # JUST_INSTALL_CEREAL: skip boost, tests and other internal stuff.
    # CEREAL_INSTALL: control if cereal gets installed.
    # Looks weird to say "just install cereal" and do not install cereal but ...
    # that's what the CMake file says.
    option(JUST_INSTALL_CEREAL "" ON)
    option(CEREAL_INSTALL "" OFF)

    # Add cereal project.
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cereal/)


endif()
