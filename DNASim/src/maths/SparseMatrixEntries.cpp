// SparseMatrixEntries class
// Nicolas Clauvelin


#include "maths/Matrix3.h"
#include "maths/Matrix4.h"
#include "maths/MatrixN.h"
#include "maths/SparseMatrixEntries.h"


namespace DNASim {


    // class constructor with size initialization
    SparseMatrixEntries::SparseMatrixEntries(Size n_entries) :
    m_entries() {
        m_entries.reserve(n_entries);
    };


    // single entry adding method
    void SparseMatrixEntries::add_entry(const EigenTriplet& entry) {
        m_entries.push_back(entry);
    };


    // entries from matrix adding method (3x3)
    void SparseMatrixEntries::add_entries_from_matrix(const Matrix3& mat3,
                                                      const Size& row_start,
                                                      const Size& col_start) {
#define MAT_SIZE 3
        // we iterate first over the columns (Eigen is column storage)
        for (Size i_col=0; i_col<MAT_SIZE; ++i_col)
            for (Size j_row=0; j_row<MAT_SIZE; ++j_row)
                m_entries.push_back(EigenTriplet((Integer)j_row+row_start,
                                                 (Integer)i_col+col_start,
                                                 mat3(j_row,i_col)));
#undef MAT_SIZE
    };


    // entries from matrix adding method (4x4)
    void SparseMatrixEntries::add_entries_from_matrix(const Matrix4& mat4,
                                                      const Size& row_start,
                                                      const Size& col_start) {
#define MAT_SIZE 4
        // we iterate first over the columns (Eigen is column storage)
        for (Size i_col=0; i_col<MAT_SIZE; ++i_col)
            for (Size j_row=0; j_row<MAT_SIZE; ++j_row)
                m_entries.push_back(EigenTriplet((Integer)j_row+row_start,
                                                 (Integer)i_col+col_start,
                                                 mat4(j_row,i_col)));
#undef MAT_SIZE
    };


    // entries from matrix adding method (nxn)
    void SparseMatrixEntries::add_entries_from_matrix(const MatrixN& mat,
                                                      const Size& row_start,
                                                      const Size& col_start) {
        const Size n_rows = mat.n_rows();
        const Size n_cols = mat.n_cols();
        for (Size i_col=0; i_col<n_cols; ++i_col)
            for (Size j_row=0; j_row<n_rows; ++j_row)
                m_entries.push_back(EigenTriplet((Integer)j_row+row_start,
                                                 (Integer)i_col+col_start,
                                                 mat(j_row,i_col)));
    };


    // clear method
    void SparseMatrixEntries::clear() {
        m_entries.clear();
    };


    // iterators
    std::vector<EigenTriplet>::const_iterator SparseMatrixEntries::begin() const {
        return m_entries.begin();
    };
    std::vector<EigenTriplet>::const_iterator SparseMatrixEntries::end() const {
        return m_entries.end();
    };


}
