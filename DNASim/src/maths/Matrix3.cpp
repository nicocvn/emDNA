// Matrix3 class
// Nicolas Clauvelin


#include "maths/Vector3.h"
#include "maths/Matrix3.h"


namespace DNASim {

    // class default constructor
    Matrix3::Matrix3() : m_mat(EigenMatrix3::Zero()) {};


    // class constructor by moving
    Matrix3::Matrix3(Matrix3&& m3) : Matrix3() {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m3.m_mat);
    };


    // instanciation methods
    Matrix3 Matrix3::random_matrix() {
        return Matrix3(EigenMatrix3::Random());
    };
    Matrix3 Matrix3::identity_matrix() {
        return Matrix3(EigenMatrix3::Identity());
    };
    Matrix3 Matrix3::skew_matrix_from_vector(const Vector3& v) {
        Matrix3 m;
        m.m_mat << FLOAT_INIT, -v[Z], v[Y], // comma initialiazer is row-order
            v[Z], FLOAT_INIT, -v[X],
            -v[Y], v[X], FLOAT_INIT;
        return m;
    };

    
    // matrix entries operators
    const Real& Matrix3::operator()(Size i, Size j) const {
        return m_mat((Integer)i,(Integer)j);
    };
    Real& Matrix3::operator()(Size i, Size j) {
        return m_mat((Integer)i,(Integer)j);
    };


    // matrix row accessor/modifier
    Vector3 Matrix3::row(Size i) const {
        EigenVector3 row_vec(m_mat.row((Integer)i));
        return Vector3(row_vec.x(), row_vec.y(), row_vec.z());
    };
    void Matrix3::set_row(Size i, const Vector3& row_vector) {
        m_mat.row((Integer)i) =
        EigenVector3(row_vector[X], row_vector[Y], row_vector[Z]);
    };

    
    // matrix col accessor/modifier
    Vector3 Matrix3::col(Size i) const {
        EigenVector3 col_vec(m_mat.col((Integer)i));
        return Vector3(col_vec.x(), col_vec.y(), col_vec.z());
    };
    void Matrix3::set_col(Size i, const Vector3& col_vector) {
        m_mat.col((Integer)i) =
        EigenVector3(col_vector[X], col_vector[Y], col_vector[Z]);
    };


    // matrix manipulations
    void Matrix3::clear() {
        m_mat.setZero();
    };
    void Matrix3::transpose() {
        m_mat.transposeInPlace();
    };
    void Matrix3::inverse() {
        EigenMatrix3 m(m_mat);
        m_mat = m.inverse();
    };
    Matrix3 Matrix3::get_transpose() const {
        EigenMatrix3 m(m_mat);
        m.transposeInPlace();
        return Matrix3(m);
    };


    // matrix eigenvalues and eigenvectors
    // values and vectors can be complex depending on the matrix
    // conditioning; the results are therefore returned as pairs of real
    // numbers (first element is the real part, second is the imaginary
    // part)
    Eigenvalues Matrix3::eigenvalues() const {

        // eigensolver
        // the second argument is a boolean flag to compute eigenvectors or not
        Eigen::EigenSolver<EigenMatrix3> solver(m_mat, Eigen::EigenvaluesOnly);

        // eigenvalues
        std::vector<RealPair> eigenvalues({
            RealPair(solver.eigenvalues()[0].real(),
                     solver.eigenvalues()[0].imag()),
            RealPair(solver.eigenvalues()[1].real(),
                     solver.eigenvalues()[1].imag()),
            RealPair(solver.eigenvalues()[2].real(),
                     solver.eigenvalues()[2].imag())
        });

        // sorting
        std::sort(eigenvalues.begin(), eigenvalues.end(),
                  [](const RealPair& p1, const RealPair& p2)->bool{
                      return p1.first > p2.first;
                  });

        return eigenvalues;

    };
    Eigensystem Matrix3::eigensystem() const {

        // eigensolver
        // the second argument is a boolean flag to compute eigenvectors or not
        Eigen::EigenSolver<EigenMatrix3> solver(m_mat,
                                                Eigen::ComputeEigenvectors);

        // eigenvalues
        std::vector<RealPair> eigenvalues({
            RealPair(solver.eigenvalues()[0].real(),
                     solver.eigenvalues()[0].imag()),
            RealPair(solver.eigenvalues()[1].real(),
                     solver.eigenvalues()[1].imag()),
            RealPair(solver.eigenvalues()[2].real(),
                     solver.eigenvalues()[2].imag())
        });

        // find sorting sequence for sorting eigenvectors
        std::vector<std::pair<Real,Size>> indexed_real_parts({
            std::pair<Real,Size>(eigenvalues[0].first, 0),
            std::pair<Real,Size>(eigenvalues[1].first, 1),
            std::pair<Real,Size>(eigenvalues[2].first, 2)
        });
        std::sort(indexed_real_parts.begin(), indexed_real_parts.end(),
                  [](const std::pair<Real,Size>& p1,
                     const std::pair<Real,Size>& p2)->bool{
                      return p1.first > p2.first;
                  });
        std::vector<Size> sorting_sequence({
            indexed_real_parts[0].second,
            indexed_real_parts[1].second,
            indexed_real_parts[2].second
        });

        // sorting
        std::sort(eigenvalues.begin(), eigenvalues.end(),
                  [](const RealPair& p1, const RealPair& p2)->bool{
                      return p1.first > p2.first;
                  });

        // eigenvectors
        Eigen::Matrix3cd vectors_mat = solver.eigenvectors();

        // conversion to pairs of real numbers
        Eigenvectors eigenvectors;
        eigenvectors.reserve(3);
        for (Size i=0; i<vectors_mat.cols(); ++i) {
            eigenvectors.push_back(std::vector<RealPair>({
                RealPair(vectors_mat.col(i)[0].real(),
                         vectors_mat.col(i)[0].imag()),
                RealPair(vectors_mat.col(i)[1].real(),
                         vectors_mat.col(i)[1].imag()),
                RealPair(vectors_mat.col(i)[2].real(),
                         vectors_mat.col(i)[2].imag())
            }));
        };

        // we need to sort the eigenvectors similarly to the eigenvalues
        Eigenvectors sorted_vectors;
        for (Size i=0; i<sorting_sequence.size(); ++i)
            sorted_vectors.push_back(eigenvectors[sorting_sequence[i]]);

        return Eigensystem(eigenvalues, sorted_vectors);

    };


    // matrix-filling operator
    Matrix3& Matrix3::operator<<(const Array<Real>& array) {

        // size checking
        DS_ASSERT(array._size == 9,
                  "Matrix3::operator<< called with wrong size;"
				  " size=" << array._size);

        // filling
        // comma initialiazer is row-order
        m_mat << array._values[0], array._values[1], array._values[2],
            array._values[3], array._values[4], array._values[5],
            array._values[6], array._values[7], array._values[8];

        return *this;
        
    };


    // block methods
    EigenArrayXX Matrix3::block_matrix() const {
        return m_mat.block(0, 0, 3, 3);
    };


    // operators
    Matrix3& Matrix3::operator=(Matrix3&& m) {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m.m_mat);
        return *this;
    };
    Matrix3& Matrix3::operator+=(const Matrix3& m)  {
        m_mat += m.m_mat;
        return *this;
    };
    Matrix3& Matrix3::operator-=(const Matrix3& m) {
        m_mat -= m.m_mat;
        return *this;
    };
    Matrix3& Matrix3::operator*=(const Real& a)  {
        m_mat *= a;
        return *this;
    };
    Matrix3& Matrix3::operator*=(const Matrix3& m)  {
        m_mat *= m.m_mat;
        return *this;
    };

    
    // matrix operators
    Matrix3 operator+(const Matrix3& m1, const Matrix3& m2) {
        return Matrix3(m1.m_mat+m2.m_mat);
    };
    Matrix3 operator-(const Matrix3& m1, const Matrix3& m2) {
        return Matrix3(m1.m_mat-m2.m_mat);
    };
    Matrix3 operator*(const Real& a, const Matrix3& m) {
        return Matrix3(a*m.m_mat);
    };
    Matrix3 operator*(const Matrix3& m1, const Matrix3& m2) {
        return Matrix3(m1.m_mat*m2.m_mat);
    };


    // vector operators
    Vector3 operator*(const Matrix3& m, const Vector3& v) {
        Eigen::Map<const EigenVector3> vec_map(v.packed_values()._values);
        EigenVector3 mv(m.m_mat*vec_map);
        return Vector3(mv.x(), mv.y(), mv.z());
    };


    // output operator
    std::ostream& operator<<(std::ostream& os, const Matrix3& m) {
        os << std::fixed << std::setprecision(REAL_WIDTH);
        os << "{";
        os << "{" << m(0,0) << ", " << m(0,1) << ", " << m(0,2) << "}, ";
        os << "{" << m(1,0) << ", " << m(1,1) << ", " << m(1,2) << "}, ";
        os << "{" << m(2,0) << ", " << m(2,1) << ", " << m(2,2) << "}";
		os << "}";
		return os;
    };

    
    // private constructor with initialization from an Eigen object
    Matrix3::Matrix3(const EigenMatrix3& m3) : m_mat(m3) {};

    
}
