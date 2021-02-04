# ---------------------------------------------------------------------------- #
# emDNA tests
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# Tests sources.
set(emDNA_tests_SOURCES
    utest_AlglibGradientCheck.cpp
    utest_AlglibGradientCheckFrozen.cpp
    utest_AlglibMinimizer.cpp
    utest_BpCollection.cpp
    utest_BpCollectionElasticEnergy.cpp
    utest_BpCollectionGradient.cpp
    utest_BpCollectionHessian.cpp
    utest_BpCollectionUpdate.cpp
    utest_BpGeometryFunctions.cpp)

# Unit tests executable.
add_executable(emDNA_tests ${emDNA_tests_SOURCES})

# Linking.
target_link_libraries(emDNA_tests PUBLIC gtest_main emDNA)


# Install directives.
install(TARGETS emDNA_tests
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
