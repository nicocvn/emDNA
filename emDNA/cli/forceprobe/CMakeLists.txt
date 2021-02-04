# ---------------------------------------------------------------------------- #
# emDNA force-probe
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# Tests sources.
set(emDNA_force_probe_SOURCES
    emDNA_ForceProbe.cpp
    ForceRampMinimizer.cpp)

# Unit tests executable.
add_executable(emDNA-force-probe ${emDNA_force_probe_SOURCES})
target_include_directories(emDNA-force-probe PUBLIC .)

# Linking.
target_link_libraries(emDNA-force-probe PUBLIC emDNA)


# Install directives.
install(TARGETS emDNA-force-probe
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
