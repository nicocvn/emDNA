# ---------------------------------------------------------------------------- #
# emDNA cli app
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# Tests sources.
set(emDNA_cli_app_SOURCES
    BpCollectionInspector.cpp
    emDNA_Application.cpp
    emDNA_CommandLine.cpp
    emDNA_LogFile.cpp)

# Unit tests executable.
add_executable(emDNA-cli ${emDNA_cli_app_SOURCES})
target_include_directories(emDNA-cli PUBLIC .)

# Linking.
target_link_libraries(emDNA-cli PUBLIC emDNA tclap)


# Install directives.
install(TARGETS emDNA-cli
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
