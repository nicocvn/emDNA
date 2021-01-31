# ---------------------------------------------------------------------------- #
# DNASim / src / maths
#
# Nicolas Clauvelin (n.clauvelin@gmail.com)
#
# ---------------------------------------------------------------------------- #


set(maths_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/AffineTransformation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Matrix3.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Matrix4.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MatrixN.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SparseMatrixEntries.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SparseMatrixN.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Triad.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Vector3.cpp
    ${CMAKE_CURRENT_LIST_DIR}/VectorN.cpp
    PARENT_SCOPE)
