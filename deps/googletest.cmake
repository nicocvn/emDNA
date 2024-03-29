# -------------------------------------------------------------------------- #
# Googletest dependency CMake file
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# -------------------------------------------------------------------------- #


if(NOT TARGET gtest)

    # Add Googletest project.
    set(BUILD_GMOCK OFF CACHE STRING "" FORCE)
    set(INSTALL_GTEST OFF CACHE STRING "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/googletest/)

    # Disable LTO (not supported due to old CMake).
    # Might need to be updated at some point.
    string(TOUPPER ${CMAKE_BUILD_TYPE} BuildType)
    set_property(TARGET gtest
                 PROPERTY INTERPROCEDURAL_OPTIMIZATION_${BuildType} FALSE)
    set_property(TARGET gtest_main
                 PROPERTY INTERPROCEDURAL_OPTIMIZATION_${BuildType} FALSE)

    # Prevent warnings.
    target_compile_options(gtest PRIVATE "-w")
    target_compile_options(gtest_main PRIVATE "-w")

endif()
