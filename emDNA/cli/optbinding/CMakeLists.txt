# ---------------------------------------------------------------------------- #
# emDNA pro bind
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


add_executable(emDNA-pro-bind
               emProBind.cpp
               ProteinBindingRamp.cpp
               RampMinimizer.cpp
               RampMinimizerEEDR.cpp
               RampMinimizerPulling.cpp)
target_include_directories(emDNA-pro-bind PUBLIC .)
target_link_libraries(emDNA-pro-bind PUBLIC emDNA)
install(TARGETS emDNA-pro-bind
        RUNTIME DESTINATION "bin/${CMAKE_BUILD_TYPE}")
