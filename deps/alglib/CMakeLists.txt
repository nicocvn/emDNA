# ---------------------------------------------------------------------------- #
# emDNA / src / alglib
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


set(alglib_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/alglibinternal.cpp
    ${CMAKE_CURRENT_LIST_DIR}/alglibmisc.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ap.cpp
    ${CMAKE_CURRENT_LIST_DIR}/integration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/interpolation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/linalg.cpp
    ${CMAKE_CURRENT_LIST_DIR}/optimization.cpp
    ${CMAKE_CURRENT_LIST_DIR}/solvers.cpp
    ${CMAKE_CURRENT_LIST_DIR}/specialfunctions.cpp)


# alglib static library.
add_library(alglib STATIC ${alglib_SOURCES})
target_include_directories(alglib PUBLIC ${CMAKE_CURRENT_LIST_DIR})