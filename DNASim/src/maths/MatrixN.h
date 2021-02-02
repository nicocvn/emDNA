// MatrixN class
// Nicolas Clauvelin


// m x n matrix implementation based on Eigen


#ifndef DNASim_MatrixN_h
#define DNASim_MatrixN_h


#include "Eigen_Includes.h"


namespace DNASim {


    class VectorN;


    class MatrixN {


    public:

        // constructors
        // default constructor creates an empty (0x0) matrix
        // by default all entries are set to zero
        MatrixN();
        MatrixN(const MatrixN& m) = default;
        MatrixN(MatrixN&& m);
        MatrixN(Size m, Size n);
        MatrixN(Size n);
        ~MatrixN() = default;

        // instantiation methods
        static MatrixN identity_matrix(Size n);
        static MatrixN diagonal_matrix(const VectorN& diagonal);
        static MatrixN random_matrix(Size m, Size n);

        // size accessors
        Size n_rows() const;
        Size n_cols() const;

        // matrix entries operators
        const Real& operator()(Size i, Size j) const;
        Real& operator()(Size i, Size j);

        // matrix manipulations
        // all entries set to zero upon resizing
        void clear();
        void transpose();
        void resize(Size nRows, Size nCols);
        MatrixN get_transpose() const;

        // static method for eigenvalues and eigenvectors
        // FOR SELF ADJOINT MATRIX ONLY !!!
        // for large matrices ... this can be hazardous (but should work)
        static std::vector<Real> self_adjoint_eigenvalues(const MatrixN& m);
        static void self_adjoint_eigensystem(const MatrixN& m,
                                             std::vector<Real>& eigenvalues,
                                             MatrixN& eigenvectors,
                                             const bool& sorted);

        // matrix-filling operator
        // array has to be row-order
        MatrixN operator<<(const Array<Real>& array);

        // matrix row accessor/modifier
        VectorN row(Size i) const;
        void set_row(Size i, const VectorN& row_vector);

        // matrix col accessor/modifier
        VectorN col(Size i) const;
        void set_col(Size i, const VectorN& col_vector);

        // block methods
        // the block is return as a row-ordered array
        EigenArrayXX block_matrix() const;
        EigenArrayXX block(Size row_index, Size col_index,
                           Size row_size, Size col_size) const;
        void set_block(Size row_index, Size col_index,
                       const EigenArrayXX& block);

        // operators
        MatrixN& operator=(const MatrixN& m) = default;
        MatrixN& operator=(MatrixN&& m);
        MatrixN& operator+=(const MatrixN& m);
        MatrixN& operator-=(const MatrixN& m);
        MatrixN& operator*=(Real a);
        MatrixN& operator*=(const MatrixN& m);

        // matrix operators
        friend MatrixN operator+(const MatrixN& m1, const MatrixN& m2);
        friend MatrixN operator-(const MatrixN& m1, const MatrixN& m2);
        friend MatrixN operator*(Real a, const MatrixN& m);
        friend MatrixN operator*(const MatrixN& m1, const MatrixN& m2);

        // vector operators
        friend VectorN operator*(const MatrixN& m, const VectorN& v);

        // output operator
        friend std::ostream& operator<<(std::ostream& os, const MatrixN& m);


    private:

        // private constructor with initialization from an Eigen object
        MatrixN(const EigenMatrixX& m);

        EigenMatrixX m_mat;


    };


}


#endif  // DNASim_MatrixN_h
