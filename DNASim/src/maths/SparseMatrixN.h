// SparseMatrixN class
// Nicolas Clauvelin


// m x n sparse matrix implementation based on Eigen

// this class does not provide data accessor; it should only be used for
// matrix operations


#ifndef DNASim_SparseMatrixN_h
#define DNASim_SparseMatrixN_h


#include "maths/Eigen_Includes.h"


namespace DNASim {


    class SparseMatrixEntries;
    class VectorN;
    class MatrixN;


    class SparseMatrixN {


    public:

        // constructors
        SparseMatrixN();
        SparseMatrixN(Size n_rows, Size n_cols);
        SparseMatrixN(const SparseMatrixN& sparse_matrix);
        SparseMatrixN(SparseMatrixN&& sparse_matrix);
        ~SparseMatrixN();

        // copy operator
        SparseMatrixN& operator=(const SparseMatrixN& sparse_matrix);

        // size accessors
        Size n_rows() const;
        Size n_cols() const;

        // resizing methods
        void resize(Size n_rows, Size n_cols);

        // extract method
        // create a new sparse matrix corresponding to the
        SparseMatrixN extract(Size row_start, Size col_start,
                              Size n_rows, Size n_cols) const;

        // filling method
        // the matrix needs to be correctly sized
        void set_nnz_entries(const SparseMatrixEntries& triplets);
        void set_nnz_entries_using_insert(const std::vector<Size>& col_nnz,
                                          const SparseMatrixEntries& triplets);

        // operators
        SparseMatrixN& operator=(SparseMatrixN&& m);
        SparseMatrixN& operator+=(const SparseMatrixN& m);
        SparseMatrixN& operator-=(const SparseMatrixN& m);
        SparseMatrixN& operator*=(Real a);
        SparseMatrixN& operator*=(const SparseMatrixN& m);

        // matrix operators
        friend MatrixN operator*(const SparseMatrixN& sm, const MatrixN& dm);
        friend MatrixN operator*(const MatrixN& dm, const SparseMatrixN& sm);

        // vector operators
        friend VectorN operator*(const SparseMatrixN& sm, const VectorN& vec);


    private:

        EigenSparseMatrix m_mat;


    };


}


#endif  // DNASim_SparseMatrixN_h
