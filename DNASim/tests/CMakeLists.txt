# ---------------------------------------------------------------------------- #
# DNASim tests
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# Tests sources.
set(DNASim_tests_SOURCES
    test_AffineTransformation.cpp
    #test_CollisionDetection.cpp
    test_CurveTopology.cpp
    #test_DNASim.cpp
    test_DNASimArchive.cpp
    test_DiscreteRibbon.cpp
    test_EnhancedString.cpp
    test_MatrixN.cpp
    test_PRN.cpp
    test_PointCloudKdTree.cpp
    test_SparseIndex.cpp
    test_StepArray.cpp
    test_StepParameters.cpp
    test_Triad.cpp
    test_Vector3.cpp
    test_VectorN.cpp
    )

# Unit tests executable.
add_executable(DNASim_tests ${DNASim_tests_SOURCES})

# Linking.
target_link_libraries(DNASim_tests PUBLIC gtest_main DNASim)
# todo: remove once Eigen is better managed as a dependency
target_include_directories(DNASim_tests PUBLIC ../deps)


# Install directives.
install(TARGETS DNASim_tests
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
