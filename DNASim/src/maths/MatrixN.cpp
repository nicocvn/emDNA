// MatrixN class
// Nicolas Clauvelin


#include <VectorN.h>
#include <MatrixN.h>


namespace DNASim {


    // class default constructor
    MatrixN::MatrixN() : m_mat() {};


    // class constructor by moving
    MatrixN::MatrixN(MatrixN&& m) : MatrixN() {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m.m_mat);
    };


    // class constructor with initialization as a rectangular matrix
    MatrixN::MatrixN(Size m, Size n) :
    m_mat(EigenMatrixX::Zero((Integer)m, (Integer)n)) {};


    // class constructor with initialization as a square matrix
    MatrixN::MatrixN(Size n) :
    m_mat(EigenMatrixX::Zero((Integer)n, (Integer)n)) {};


    // instantiation methods
    MatrixN MatrixN::identity_matrix(Size n) {
        return MatrixN(EigenMatrixX::Identity((Integer)n,(Integer)n));
    };
    MatrixN MatrixN::diagonal_matrix(const VectorN& diagonal) {
        Eigen::Map<const EigenVectorX>
            vec_map(diagonal.pack_values()._values, (Integer)diagonal.size());
        return MatrixN(vec_map.asDiagonal());
    };
    MatrixN MatrixN::random_matrix(Size m, Size n) {
        return MatrixN(EigenMatrixX::Random((Integer)m,(Integer)n));
    };


    // size accessors
    Size MatrixN::n_rows() const { return (Size)m_mat.rows(); };
    Size MatrixN::n_cols() const { return (Size)m_mat.cols(); };

    
    // matrix entries operators
    const Real& MatrixN::operator()(Size i, Size j) const {
        return m_mat((Integer)i,(Integer)j);
    };
    Real& MatrixN::operator()(Size i, Size j) {
        return m_mat((Integer)i,(Integer)j);
    };


    // matrix manipulations
    void MatrixN::clear() { m_mat.setZero(); };
    void MatrixN::transpose() { m_mat.transposeInPlace(); };
    void MatrixN::resize(Size nRows, Size nCols) {
        m_mat.resize((Integer)nRows, (Integer)nCols);
        m_mat.setZero();
    };
    MatrixN MatrixN::get_transpose() const {
        EigenMatrixX m(m_mat);
        m.transposeInPlace();
        return MatrixN(m);
    };


    // static method for eigenvalues and eigenvectors
    // FOR SELF ADJOINT MATRIX ONLY !!!
    // for large matrices ... this can be hazardous (but should work)
    std::vector<Real> MatrixN::self_adjoint_eigenvalues(const MatrixN& m) {

        // check if the matrix is square
        DS_ASSERT(m.n_rows()==m.n_cols(),
                  "computing eigendecomposition with non-square matrix");

        // sizing
        const Size mat_size = m.n_rows();

        // solver construction
        Eigen::SelfAdjointEigenSolver<EigenMatrixX> solver(mat_size);

        // eigenvalues computation
        // the second argument is to prevent eigenvectors computation
        solver.compute(m.m_mat, Eigen::EigenvaluesOnly);

        // eigenvalues
        std::vector<Real> eigenvalues;
        eigenvalues.reserve(mat_size);
        for (Size i=0; i<mat_size; ++i)
            eigenvalues.push_back(solver.eigenvalues().operator()(i,0));

        // sorting
        std::sort(eigenvalues.begin(), eigenvalues.end(),
                  [](const Real& x1, const Real& x2) {
                      return x1<x2;
                  });

        return eigenvalues;

    };
    void MatrixN::self_adjoint_eigensystem(const MatrixN& m,
                                           std::vector<Real>& eigenvalues,
                                           MatrixN& eigenvectors,
                                           const bool& sorted) {

        // check if the matrix is square
        DS_ASSERT(m.n_rows()==m.n_cols(),
                  "computing eigendecomposition with non-square matrix");

        // sizing
        const Size mat_size = m.n_rows();

        // cleaning and allocation
        eigenvalues.clear();
        eigenvalues.reserve(mat_size);
        eigenvectors = MatrixN(mat_size);

        // solver construction
        Eigen::SelfAdjointEigenSolver<EigenMatrixX> solver(mat_size);

        // eigensystem computation
        solver.compute(m.m_mat, Eigen::ComputeEigenvectors);

        // building data structure for sorting
        EigenMatrixX values = solver.eigenvalues();
        EigenMatrixX vectors = solver.eigenvectors();
        std::vector<EigenPair> eigen_pairs;
        eigen_pairs.reserve(mat_size);
        for (Size i=0; i<mat_size; ++i) {

            // prepare stl vector
            std::vector<Real> v(mat_size);
            for (Size k=0; k<mat_size; ++k)
                v[k] = vectors(k,i);

            eigen_pairs.push_back(EigenPair(values(i,0), v));

        };

        // sorting
        if (sorted)
            std::sort(eigen_pairs.begin(), eigen_pairs.end(),
                      [](const EigenPair& p1, const EigenPair& p2) {
                          return p1.first < p2.first;
                      });

        // unpacking
        for (Size i=0; i<mat_size; ++i) {

            // value
            eigenvalues.push_back(eigen_pairs[i].first);

            // vector
            eigenvectors.set_col(i, VectorN(eigen_pairs[i].second));

        };

    };


    // matrix-filling operator
    MatrixN MatrixN::operator<<(const Array<Real>& array) {

        // size checking
        DS_ASSERT(array._size == (Size)(m_mat.cols()*m_mat.rows()),
                  "MatrixN::operator<< called with wrong size;"
				  " size=" << array._size <<
                  " ; matrix size=" << m_mat.cols()*m_mat.rows());

        // filling
        // we have to transpose the matrix because of the internal storage
        // in Eigen object (column-order by default)
        Eigen::Map<const EigenMatrixX> mat_map(array._values,
                                               m_mat.rows(), m_mat.cols());
        m_mat.noalias() = mat_map.transpose();

        return *this;

    };


    // matrix row accessor/modifier
    VectorN MatrixN::row(Size i) const {
        EigenVectorX row_vec(m_mat.row((Integer)i));
        return VectorN(Array<Real>((Size)row_vec.size(), &row_vec(0)));
    };
    void MatrixN::set_row(Size i, const VectorN& row_vector) {
        Eigen::Map<const EigenVectorX> vec_map(row_vector.pack_values()._values,
                                               (Integer)row_vector.
                                               pack_values()._size);
        m_mat.row((Integer)i) = vec_map;
    };


    // matrix col accessor/modifier
    VectorN MatrixN::col(Size i) const {
        EigenVectorX col_vec(m_mat.col((Integer)i));
        return VectorN(Array<Real>((Size)col_vec.size(), &col_vec(0)));
    };
    void MatrixN::set_col(Size i, const VectorN& col_vector) {
        Eigen::Map<const EigenVectorX> vec_map(col_vector.pack_values()._values,
                                               (Integer)col_vector.
                                               pack_values()._size);
        m_mat.col((Integer)i) = vec_map;
    };


    // block methods
    EigenArrayXX MatrixN::block_matrix() const {
        return m_mat.block(0, 0, m_mat.rows(), m_mat.cols());
    };
    EigenArrayXX MatrixN::block(Size row_index, Size col_index,
                               Size row_size, Size col_size) const {
        return m_mat.block(row_index, col_index, row_size, col_size);
    };
    void MatrixN::set_block(Size row_index, Size col_index,
                            const EigenArrayXX& block) {
        m_mat.block(row_index, col_index,
                    block.rows(), block.cols()) = block;
    };


    // operators
    MatrixN& MatrixN::operator=(MatrixN&& m) {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m.m_mat);
        return *this;
    };
    MatrixN& MatrixN::operator+=(const MatrixN& m) {
        m_mat.noalias() += m.m_mat;
        return *this;
    };
    MatrixN& MatrixN::operator-=(const MatrixN& m) {
        m_mat.noalias() -= m.m_mat;
        return *this;
    };
    MatrixN& MatrixN::operator*=(Real a) {
        m_mat *= a;
        return *this;
    };
    MatrixN& MatrixN::operator*=(const MatrixN& m) {
        m_mat *= m.m_mat;
        return *this;
    };

    
    // matrix operators
    MatrixN operator+(const MatrixN& m1, const MatrixN& m2) {
        return MatrixN(m1.m_mat+m2.m_mat);
    };
    MatrixN operator-(const MatrixN& m1, const MatrixN& m2) {
        return MatrixN(m1.m_mat-m2.m_mat);
    };
    MatrixN operator*(Real a, const MatrixN& m) {
        return MatrixN(a*m.m_mat);
    };
    MatrixN operator*(const MatrixN& m1, const MatrixN& m2) {
        return MatrixN(m1.m_mat*m2.m_mat);
    };


    // vector operators
    VectorN operator*(const MatrixN& m, const VectorN& v) {
        Eigen::Map<const EigenVectorX> vec_map(v.pack_values()._values,
                                               (Integer)v.pack_values()._size);
        EigenVectorX mv(m.m_mat*vec_map);
        return VectorN(Array<Real>((Size)mv.size(), &(mv(0))));
    };


    // output operator
    std::ostream& operator<<(std::ostream& os, const MatrixN& m) {
        os << std::fixed << std::setprecision(REAL_WIDTH);
        Size n_rows = m.n_rows();
        Size n_cols = m.n_cols();
        os << "{";
        for (Size i=0; i<n_rows-1; ++i) {
            os << "{";
            for (Size j=0; j<n_cols-1; ++j)
                os << m(i,j) << ", ";
            os << m(i,n_cols-1) << "}, ";
        };
        os << "{";
        for (Size j=0; j<n_cols-1; ++j)
            os << m(n_rows-1,j) << ", ";
        os << m(n_rows-1,n_cols-1) << "}";
        os << "}";
        return os;
    };


    // private constructor with initialization from an Eigen object
    MatrixN::MatrixN(const EigenMatrixX& m) : m_mat(m) {};


}
