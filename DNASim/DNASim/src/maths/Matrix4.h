// Matrix4 class
// Nicolas Clauvelin


// 4x4 matrix implementation based on Eigen


#ifndef DNASim_Matrix4_h
#define DNASim_Matrix4_h


#include <Eigen_Includes.h>


namespace DNASim {


    class Matrix4 {


    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        // constructors
        // default constructor creates a zero matrix
        Matrix4();
        Matrix4(const Matrix4& m4) = default;
        Matrix4(Matrix4&& m4);
        ~Matrix4() = default;

        // instanciation methods
        static Matrix4 random_matrix();
        static Matrix4 identity_matrix();

        // matrix entries operators
        const Real& operator()(Size i, Size j) const;
        Real& operator()(Size i, Size j);
        
        // matrix manipulations
		void clear();
		void transpose();
        void inverse();
        Matrix4 get_transpose() const;

        // matrix-filling operator
        // array has to be row-order
        Matrix4& operator<<(const Array<Real>& array);

        // block method
        // return the matrix as a row-ordered block
        EigenArrayXX block_matrix() const;

        // operators
        Matrix4& operator=(const Matrix4& m) = default;
        Matrix4& operator=(Matrix4&& m);
		Matrix4& operator+=(const Matrix4& m);
		Matrix4& operator-=(const Matrix4& m);
        Matrix4& operator*=(const Real& a);
		Matrix4& operator*=(const Matrix4& m);

		// matrix operators
		friend Matrix4 operator+(const Matrix4& m1, const Matrix4& m2);
		friend Matrix4 operator-(const Matrix4& m1, const Matrix4& m2);
		friend Matrix4 operator*(const Matrix4& m1, const Matrix4& m2);
        friend Matrix4 operator*(const Real& a, const Matrix4& m);

        // output operator
        friend std::ostream& operator<<(std::ostream& os, const Matrix4& m);


    private:

        // private constructor with initialization from an Eigen object
        explicit Matrix4(const EigenMatrix4& m4);
        
        EigenMatrix4 m_mat;
        
        
    };
    
    
}



#endif  // DNASim_Matrix4_h
