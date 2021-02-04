# ---------------------------------------------------------------------------- #
# emDNA tools
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


# check-collision
add_executable(emDNA-check-collision emDNA_check_collision.cpp)
target_include_directories(emDNA-check-collision PUBLIC .)
target_link_libraries(emDNA-check-collision PUBLIC emDNA tclap)
install(TARGETS emDNA-check-collision
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")


# lego-gen
add_executable(emDNA-lego-gen emDNA_lego_gen.cpp)
target_include_directories(emDNA-lego-gen PUBLIC .)
target_link_libraries(emDNA-lego-gen PUBLIC emDNA tclap)
install(TARGETS emDNA-lego-gen
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")


# parser
add_executable(emDNA-parser emDNA_parser.cpp)
target_include_directories(emDNA-parser PUBLIC .)
target_link_libraries(emDNA-parser PUBLIC emDNA tclap)
install(TARGETS emDNA-parser
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")


# topology
add_executable(emDNA-topology emDNA_topology.cpp)
target_include_directories(emDNA-topology PUBLIC .)
target_link_libraries(emDNA-topology PUBLIC emDNA tclap)
install(TARGETS emDNA-topology
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
