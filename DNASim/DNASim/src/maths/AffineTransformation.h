// AffineTransformation class
// Nicolas Clauvelin


// 3-dimensional affine transformation implementation based on Eigen


#ifndef DNASim_AffineTransformation_h
#define DNASim_AffineTransformation_h


#include <Eigen_Includes.h>


namespace DNASim {


    class Matrix4;
    class Matrix3;
    class Vector3;


    class AffineTransformation {


    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        // constructors
        // default constructor creates an identity transformation
        AffineTransformation() = default;
        AffineTransformation(const AffineTransformation& aff_t) = default;
        AffineTransformation(AffineTransformation&& aff_t);
        AffineTransformation(const Matrix4& mat);
        ~AffineTransformation() = default;

        // copy and move operators
        AffineTransformation&
        operator=(const AffineTransformation& aff_t) = default;
        AffineTransformation&
        operator=(AffineTransformation&& aff_t);

        // instanciaton methods
        static AffineTransformation rotation(const Real& angle,
                                             const Vector3& axis);
        static AffineTransformation translation(const Vector3& t_vec);

        // transformation components accessors
        Matrix4 matrix() const;
        Matrix3 rotation() const;
        Vector3 translation() const;

        // transformation concatenation methods
        void append(const AffineTransformation& aff_t);
        void prepend(const AffineTransformation& aff_t);

        // inverse method
        AffineTransformation inverse() const;

        // transformation method
        Vector3 transform_point(const Vector3& point) const;
        Vector3 transform_vector(const Vector3& vector) const;


    private:

        // private constructor with initialization from an Eigen object
        AffineTransformation(const EigenTransform3D& affine_trans);

        EigenTransform3D m_trans;


    };


}



#endif  // DNASim_AffineTransformation_h
