# ---------------------------------------------------------------------------- #
# emDNA force-field-packager
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


add_executable(emDNA-ff-packager ffield_packager.cpp)
target_include_directories(emDNA-ff-packager PUBLIC .)
target_link_libraries(emDNA-ff-packager PUBLIC emDNA tclap)
install(TARGETS emDNA-ff-packager
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
