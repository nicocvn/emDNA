# ---------------------------------------------------------------------------- #
# DNASim / src / dna
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


set(dna_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/BpGeometry.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ForceConstantsDB.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Parser_x3DNA.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Sequence.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SequenceDepenceModels.cpp
    ${CMAKE_CURRENT_LIST_DIR}/StepArray.cpp
    ${CMAKE_CURRENT_LIST_DIR}/StepParameters.cpp
    ${CMAKE_CURRENT_LIST_DIR}/StepParametersDB.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TetramerDepenceModels.cpp
    PARENT_SCOPE)
