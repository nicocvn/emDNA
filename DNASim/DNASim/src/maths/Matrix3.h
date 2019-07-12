// Matrix3 class
// Nicolas Clauvelin


// 3x3 matrix implementation based on Eigen


#ifndef DNASim_Matrix3_h
#define DNASim_Matrix3_h


#include <Eigen_Includes.h>


namespace DNASim {


    class Vector3;


    class Matrix3 {


    public:

        // constructors
        // default constructor creates a zero matrix
        Matrix3();
        Matrix3(const Matrix3& m3) = default;
        Matrix3(Matrix3&& m3);
        ~Matrix3() = default;

        // instanciation methods
        static Matrix3 random_matrix();
        static Matrix3 identity_matrix();
        static Matrix3 skew_matrix_from_vector(const Vector3& v);

        // matrix entries operators
        const Real& operator()(Size i, Size j) const;
        Real& operator()(Size i, Size j);

        // matrix row accessor/modifier
        Vector3 row(Size i) const;
        void set_row(Size i, const Vector3& row_vector);

        // matrix col accessor/modifier
        Vector3 col(Size i) const;
        void set_col(Size i, const Vector3& col_vector);

        // matrix manipulations
		void clear();
		void transpose();
        void inverse();
        Matrix3 get_transpose() const;

        // matrix eigenvalues and eigenvectors
        // values and vectors can be complex depending on the matrix
        // conditioning; the results are therefore returned as pairs of real
        // numbers (first element is the real part, second is the imaginary
        // part)
        Eigenvalues eigenvalues() const;
        Eigensystem eigensystem() const;

        // matrix-filling operator
        // array has to be row-order
        Matrix3& operator<<(const Array<Real>& array);

        // block methods
        // return the matrix as a row-ordered block
        EigenArrayXX block_matrix() const;

        // operators
        Matrix3& operator=(const Matrix3& m) = default;
        Matrix3& operator=(Matrix3&& m);
		Matrix3& operator+=(const Matrix3& m);
		Matrix3& operator-=(const Matrix3& m);
        Matrix3& operator*=(const Real& a);
		Matrix3& operator*=(const Matrix3& m);

		// matrix operators
		friend Matrix3 operator+(const Matrix3& m1, const Matrix3& m2);
		friend Matrix3 operator-(const Matrix3& m1, const Matrix3& m2);
		friend Matrix3 operator*(const Matrix3& m1, const Matrix3& m2);
        friend Matrix3 operator*(const Real& a, const Matrix3& m);

        // vector operators
		friend Vector3 operator*(const Matrix3& m, const Vector3& v);

        // output operator
        friend std::ostream& operator<<(std::ostream& os, const Matrix3& m);
        

    private:

        // private constructor with initialization from an Eigen object
        explicit Matrix3(const EigenMatrix3& m3);

        EigenMatrix3 m_mat;


    };


}



#endif  // DNASim_Matrix3_h
