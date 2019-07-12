// AffineTransformation class
// Nicolas Clauvelin


#include <Matrix4.h>
#include <Matrix3.h>
#include <Vector3.h>
#include <AffineTransformation.h>


namespace DNASim {


    // class constructor with initialization from a 4x4 matrix
    AffineTransformation::AffineTransformation(const Matrix4& mat) : m_trans() {
        EigenMatrix4 m;
        m << mat(0,0), mat(0,1), mat(0,2), mat(0,3), // comma initializer
            mat(1,0), mat(1,1), mat(1,2), mat(1,3),  // is row-order
            mat(2,0), mat(2,1), mat(2,2), mat(2,3),
            mat(3,0), mat(3,1), mat(3,2), mat(3,3);
        m_trans.matrix() = m;
    };


    // constructor by moving
    AffineTransformation::AffineTransformation(AffineTransformation&& aff_t) :
    m_trans() {
        // Eigen types are not movable at this point so we rely on the swap
        m_trans.matrix().swap(aff_t.m_trans.matrix());
    };


    // move operator
    AffineTransformation& AffineTransformation::operator=(AffineTransformation&&
                                                          aff_t) {
        // Eigen types are not movable at this point so we rely on the swap
        m_trans.matrix().swap(aff_t.m_trans.matrix());
        return *this;
    };


    // instanciation methods
    AffineTransformation AffineTransformation::rotation(const Real& angle,
                                                        const Vector3& axis) {
        // normalized axis
        EigenVector3 v(axis[X], axis[Y], axis[Z]);
        v.normalize();

        // rotation transformation as quaternion
        // Eigen doc: for a single vector transformation quaternion is the
        // preferred representation
        AffineTransformation t;
        t.m_trans = Eigen::AngleAxis<Real>(angle, v);

        return t;

    };
    AffineTransformation AffineTransformation::translation(const Vector3&
                                                           t_vec) {
        Eigen::Map<const EigenVector3> vec_map(&(t_vec.packed_values().
                                                 _values[0]));
        AffineTransformation t;
        t.m_trans = Eigen::Translation<Real,3>(vec_map);
        return t;
    };
    

    // transformation components accessors
    Matrix4 AffineTransformation::matrix() const {

        // we use a tranpose because the Array filling method is row-order
        const EigenMatrix4 mat(m_trans.matrix().transpose());
        Matrix4 m;
        m << Array<Real>(16, mat.data());

        return m;
    };
    Matrix3 AffineTransformation::rotation() const {

        // we use a tranpose because the Array filling method is row-order
        const EigenMatrix3 mat(m_trans.linear().transpose());
        Matrix3 m;
        m << Array<Real>(9, mat.data());

        return m;

    };
    Vector3 AffineTransformation::translation() const {
        const EigenVector3 v(m_trans.translation());
        return Vector3(v[0], v[1], v[2]);
    };


    // transformation concatenation methods
    void AffineTransformation::append(const AffineTransformation& aff_t) {
        m_trans = m_trans*aff_t.m_trans;
    };
    void AffineTransformation::prepend(const AffineTransformation& aff_t) {
        m_trans = aff_t.m_trans*m_trans;
    };


    // inverse method
    AffineTransformation AffineTransformation::inverse() const {
        return AffineTransformation(m_trans.inverse());
    };


    // transformation method
    Vector3 AffineTransformation::transform_point(const Vector3& point) const {
        const Eigen::Map<const EigenVector3> vec_map(&(point.packed_values().
                                                       _values[0]));
        const EigenVector3 tp = m_trans*vec_map;
        return Vector3(tp[0], tp[1], tp[2]);
    };
    Vector3 AffineTransformation::transform_vector(const Vector3& vector) const
    {
        const Eigen::Map<const EigenVector3> vec_map(&(vector.packed_values().
                                                       _values[0]));
        const EigenVector3 tp = (m_trans.linear())*vec_map;
        return Vector3(tp[0], tp[1], tp[2]);
    };


    // private constructor with initialization from an Eigen object
    AffineTransformation::AffineTransformation(const EigenTransform3D&
                                               affine_trans) :
        m_trans(affine_trans) {};


}
