// SparseMatrixN class
// Nicolas Clauvelin


#include <SparseMatrixEntries.h>
#include <MatrixN.h>
#include <VectorN.h>
#include <SparseMatrixN.h>


namespace DNASim {


    // class default constructor
    SparseMatrixN::SparseMatrixN() : m_mat() {};


    // class constructor with size initialization
    SparseMatrixN::SparseMatrixN(Size n_rows, Size n_cols) :
    m_mat(n_rows, n_cols) {};


    // class constructor by copy
    SparseMatrixN::SparseMatrixN(const SparseMatrixN& sparse_matrix) :
    m_mat(sparse_matrix.m_mat) {};


    // class constructor by moving
    SparseMatrixN::SparseMatrixN(SparseMatrixN&& sparse_matrix) :
    SparseMatrixN() {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(sparse_matrix.m_mat);
    };


    // class destructor
    SparseMatrixN::~SparseMatrixN() {};


    // copy operator
    SparseMatrixN& SparseMatrixN::operator=(const SparseMatrixN& sparse_matrix)
    {
        m_mat = sparse_matrix.m_mat;
        return *this;
    };


    // size accessors
    Size SparseMatrixN::n_rows() const {
        return m_mat.rows();
    };
    Size SparseMatrixN::n_cols() const {
        return m_mat.cols();
    };


    // resizing methods
    void SparseMatrixN::resize(Size n_rows, Size n_cols) {
        m_mat.resize(n_rows, n_cols);
    };


    // extract method
    // create a new sparse matrix corresponding to the
    SparseMatrixN SparseMatrixN::extract(Size row_start, Size col_start,
                                         Size n_rows, Size n_cols) const {

        // sparse matrix entries
        SparseMatrixEntries entries(m_mat.nonZeros());

        // iterate over non-zero elements
        for (Integer k=0; k<m_mat.outerSize(); ++k)
            for (EigenSparseMatrix::InnerIterator it(m_mat,k); it; ++it) {

                // check the current non-zero element belongs to the block
                bool row_check((it.row() < (Integer)(row_start+n_rows)) &&
                               (it.row() >= (Integer)row_start));
                bool col_check((it.col() < (Integer)(col_start+n_cols)) &&
                               (it.col() >= (Integer)col_start));
                if (row_check && col_check)
                    entries.add_entry(EigenTriplet(it.row()-row_start,
                                                   it.col()-col_start,
                                                   it.value()));
            };

        // new sparse matrix
        SparseMatrixN mat(n_rows, n_cols);
        mat.set_nnz_entries(entries);

        return mat;

    };


    // filling method (from triplets)
    void SparseMatrixN::set_nnz_entries(const SparseMatrixEntries& triplets) {
        m_mat.setFromTriplets(triplets.begin(), triplets.end());
    };


    // filling method (using insert)
    void SparseMatrixN::
    set_nnz_entries_using_insert(const std::vector<Size>& col_nnz,
                                 const SparseMatrixEntries& triplets) {

        // sizing
        Eigen::VectorXi eigen_col_nnz(col_nnz.size());
        for (Size i=0, n=col_nnz.size(); i<n; ++i)
            eigen_col_nnz(i) = col_nnz[i];
        m_mat.reserve(eigen_col_nnz);

        // filling
        // optimal if the triplets is "column ordered"
        for (const EigenTriplet& triplet : triplets)
            m_mat.insert(triplet.row(), triplet.col()) = triplet.value();

        // compression
        // diabled cause we assume that the sizing is perfect !
        //m_mat.makeCompressed();

    };


    // operators
    SparseMatrixN& SparseMatrixN::operator=(SparseMatrixN&& m) {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m.m_mat);
        return *this;
    };
    SparseMatrixN& SparseMatrixN::operator+=(const SparseMatrixN& m) {
        m_mat += m.m_mat;
        return *this;
    };
    SparseMatrixN& SparseMatrixN::operator-=(const SparseMatrixN& m) {
        m_mat -= m.m_mat;
        return *this;
    };
    SparseMatrixN& SparseMatrixN::operator*=(Real a) {
        m_mat *= a;
        return *this;
    };
    SparseMatrixN& SparseMatrixN::operator*=(const SparseMatrixN& m) {
        m_mat = m_mat*m.m_mat;
        return *this;
    };


    // matrix operators
    MatrixN operator*(const SparseMatrixN& sm, const MatrixN& dm) {
        EigenMatrixX dmat(dm.block_matrix());
        MatrixN final_mat(sm.n_rows(), dm.n_cols());
        final_mat.set_block(0, 0, sm.m_mat*dmat);
        return final_mat;
    };
    MatrixN operator*(const MatrixN& dm, const SparseMatrixN& sm) {
        EigenMatrixX dmat(dm.block_matrix());
        MatrixN final_mat(dm.n_rows(), sm.n_cols());
        final_mat.set_block(0, 0, dmat*sm.m_mat);
        return final_mat;
    };


    // vector operators
    VectorN operator*(const SparseMatrixN& sm, const VectorN& vec) {
        Eigen::Map<const EigenVectorX> vec_map(vec.pack_values()._values,
                                               (Integer)vec.pack_values().
                                               _size);
        EigenVectorX mv(sm.m_mat*vec_map);
        return VectorN(Array<Real>((Size)mv.size(), &(mv(0))));
    };


}
