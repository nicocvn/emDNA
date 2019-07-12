// Matrix4 class
// Nicolas Clauvelin


#include <Matrix4.h>


namespace DNASim {

    // class default constructor
    Matrix4::Matrix4() : m_mat(EigenMatrix4::Zero()) {};


    // class constructor by moving
    Matrix4::Matrix4(Matrix4&& m4) : Matrix4() {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m4.m_mat);
    };


    // instanciation methods
    Matrix4 Matrix4::random_matrix() {
        return Matrix4(EigenMatrix4::Random());
    };
    Matrix4 Matrix4::identity_matrix() {
        return Matrix4(EigenMatrix4::Identity());
    };


    // matrix entries operators
    const Real& Matrix4::operator()(Size i, Size j) const {
        return m_mat((Integer)i,(Integer)j);
    };
    Real& Matrix4::operator()(Size i, Size j) {
        return m_mat((Integer)i,(Integer)j);
    };


    // matrix manipulations
    void Matrix4::clear() {
        m_mat.setZero();
    };
    void Matrix4::transpose() {
        m_mat.transposeInPlace();
    };
    void Matrix4::inverse() {
        EigenMatrix4 m(m_mat);
        m_mat = m.inverse();
    };
    Matrix4 Matrix4::get_transpose() const {
        EigenMatrix4 m(m_mat);
        m.transposeInPlace();
        return Matrix4(m);
    };


    // matrix-filling operator
    Matrix4& Matrix4::operator<<(const Array<Real>& array) {

        // size checking
        DS_ASSERT(array._size == 16,
                  "Matrix4::operator<< called with wrong size;"
				  " size=" << array._size);

        // filling
        // comma initialiazer is row-order
        m_mat << array._values[0], array._values[1], array._values[2],
            array._values[3], array._values[4], array._values[5],
            array._values[6], array._values[7], array._values[8],
            array._values[9], array._values[10], array._values[11],
            array._values[12], array._values[13], array._values[14],
            array._values[15];

        return *this;
        
    };


    // block methods
    EigenArrayXX Matrix4::block_matrix() const {
        return m_mat.block(0, 0, 4, 4);
    };


    // operators
    Matrix4& Matrix4::operator=(Matrix4&& m) {
        // Eigen types are not movable at this point so we rely on the swap
        m_mat.swap(m.m_mat);
        return *this;
    };
    Matrix4& Matrix4::operator+=(const Matrix4& m)  {
        m_mat += m.m_mat;
        return *this;
    };
    Matrix4& Matrix4::operator-=(const Matrix4& m) {
        m_mat -= m.m_mat;
        return *this;
    };
    Matrix4& Matrix4::operator*=(const Real& a)  {
        m_mat *= a;
        return *this;
    };
    Matrix4& Matrix4::operator*=(const Matrix4& m)  {
        m_mat *= m.m_mat;
        return *this;
    };


    // matrix operators
    Matrix4 operator+(const Matrix4& m1, const Matrix4& m2) {
        return Matrix4(m1.m_mat+m2.m_mat);
    };
    Matrix4 operator-(const Matrix4& m1, const Matrix4& m2) {
        return Matrix4(m1.m_mat-m2.m_mat);
    };
    Matrix4 operator*(const Real& a, const Matrix4& m) {
        return Matrix4(a*m.m_mat);
    };
    Matrix4 operator*(const Matrix4& m1, const Matrix4& m2) {
        return Matrix4(m1.m_mat*m2.m_mat);
    };


    // output operator
    std::ostream& operator<<(std::ostream& os, const Matrix4& m) {
        os << std::fixed << std::setprecision(REAL_WIDTH);
        os << "{";
        os << "{" << m(0,0) << ", " << m(0,1) << ", " << m(0,2)
            << ", " << m(0,3) << "}, ";
        os << "{" << m(1,0) << ", " << m(1,1) << ", " << m(1,2)
            << ", " << m(1,3) << "}, ";
        os << "{" << m(2,0) << ", " << m(2,1) << ", " << m(2,2)
            << ", " << m(2,3) << "}, ";
        os << "{" << m(3,0) << ", " << m(3,1) << ", " << m(3,2)
            << ", " << m(3,3) << "}";
		os << "}";
		return os;
    };


    // private constructor with initialization from an Eigen object
    Matrix4::Matrix4(const EigenMatrix4& m4) : m_mat(m4) {};
    
    
}
