# ---------------------------------------------------------------------------- #
# DNASim / src / pdb
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


set(pdb_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/PDBAtom.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PDBChain.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PDBParser.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PDBResidue.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PDBStructure.cpp
    PARENT_SCOPE)
