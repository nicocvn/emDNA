// Vector3 class header file
// Nicolas Clauvelin

// three-dimensional Real vector implementation based on Eigen
// the vector is implemented as a 4D vector to take advantage of vectorization


#ifndef DNASim_Vector3_h
#define DNASim_Vector3_h


#include <Eigen_Includes.h>


// Cartesian coordinates enum
enum CartesianCoordinates {
    X = 0,
    Y = 1,
    Z = 2
};


namespace DNASim {


    class Vector3 {


    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        // constructors
        Vector3();
        Vector3(const Vector3& v) = default;
        Vector3(Vector3&& v);
        Vector3(const Real& x, const Real& y, const Real& z);
        Vector3(const std::initializer_list<Real>& v);
        explicit Vector3(const std::string& s);
        ~Vector3() = default;

        // static instanciation methods
        static Vector3 random_vector();
        static Vector3 unit_x();
        static Vector3 unit_y();
        static Vector3 unit_z();

        // component const-accessor
        inline const Real& operator[](const CartesianCoordinates& c) const {
            return m_vec[static_cast<Size>(c)];
        };
        // component accessor
        inline Real& operator[](const CartesianCoordinates& c) {
            return m_vec[static_cast<Size>(c)];
        };

        // vector norm methods
        Real norm() const;
        void normalize();

        // vector operation methods
        Real dot(const Vector3& v) const;
        Vector3 cross(const Vector3& v) const;

        // vector rotation method
        // the rotation of the vector is centered at the origin of the global
        // reference frame
        // this method should not be used to rotate points with respect to a
        // specific frame
        // (see AffineTransformation class)
        void rotate(const Real& angle, const Vector3& axis);

        // packed values methods
        Array<Real> packed_values() const;
        std::vector<Real> std_vector() const;

        // operators
        Vector3& operator=(const Vector3& v) = default;
        Vector3& operator=(Vector3&& v);
        Vector3& operator+=(const Vector3& v);
        Vector3& operator-=(const Vector3& v);
        Vector3& operator*=(const Real& a);

        // vector operators
        friend Vector3 operator+(const Vector3& v1, const Vector3& v2);
        friend Vector3 operator-(const Vector3& v1, const Vector3& v2);
        friend Vector3 operator*(const Vector3& v, const Real& a);
        friend Vector3 operator*(const Real& a, const Vector3& v);

        // output operator
        friend std::ostream& operator<<(std::ostream& os, const Vector3& v);


    private:

        // private constructor with initialization from an Eigen object
        explicit Vector3(const EigenVector4& v);

        EigenVector4 m_vec;


    };


}


#endif  // DNASim_Vector3_h
