# ---------------------------------------------------------------------------- #
# emDNA / src / minim
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


set(minim_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/AlglibBpCollectionMinimizer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Alglib_Includes.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MinimizerAgent.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Minimizer_AlglibCG.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Minimizer_AlglibLBFGS.cpp
    PARENT_SCOPE)
