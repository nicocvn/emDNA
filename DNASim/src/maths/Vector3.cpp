// Vector3 class implementation file
// Nicolas Clauvelin


#include <EnhancedString.h>
#include <Vector3.h>


namespace DNASim {


    // default constructor
    Vector3::Vector3() : m_vec(FLOAT_INIT, FLOAT_INIT, FLOAT_INIT, FLOAT_INIT)
    {};


    // constructor by moving
    Vector3::Vector3(Vector3&& v) : Vector3() {
        // Eigen types are not movable at this point so we rely on the swap
        m_vec.swap(v.m_vec);
    };


    // constructor with components initialization
    Vector3::Vector3(const Real& x, const Real& y, const Real& z) :
        m_vec(x, y, z, FLOAT_INIT) {};


    // constructor with initalization from a initializer_list
    Vector3::Vector3(const std::initializer_list<Real>& v) : Vector3() {
        DS_ASSERT(v.size()==3,
                  "Vector3(std::initializer_list) called with a badly sized "
                  "list; size=" << v.size());
        m_vec[0] = *(v.begin());
        m_vec[1] = *(v.begin()+1);
        m_vec[2] = *(v.begin()+2);
    };


    // constructor with initialization from a string
    // the string is expected to be of the form {X, Y, Z}
    Vector3::Vector3(const std::string& s) : Vector3() {

        // string cleaning
		EnhancedString es(s);
		es.erase_blank_characters();
        es.erase_characters({'{','}'});

		// tokenizing
		std::vector<EnhancedString> tokens;
		es.tokenize(tokens, ',');
		DS_ASSERT(tokens.size()==3,
				  "Vector3(std::string) called with wrong size;"
				  " size=" << tokens.size());

        // filling
		m_vec[0] = tokens[0].convert<Real>();
		m_vec[1] = tokens[1].convert<Real>();
		m_vec[2] = tokens[2].convert<Real>();

    };


    // static instantiation method as a random vector
    Vector3 Vector3::random_vector() {
        return Vector3(EigenVector4::Random());
    };

    // static instantiation method as {1,0,0}
    Vector3 Vector3::unit_x() {
        return Vector3(EigenVector4::UnitX());
    };

    // static instantiation method as {0,1,0}
    Vector3 Vector3::unit_y() {
        return Vector3(EigenVector4::UnitY());
    };

    // static instantiation method as {0,0,1}
    Vector3 Vector3::unit_z() {
        return Vector3(EigenVector4::UnitZ());
    };

    
    // norm method
    Real Vector3::norm() const {
        return m_vec.norm();
    };

    // normalize method
    void Vector3::normalize() {
        m_vec.normalize();
    };


    // scalar product computation method
    Real Vector3::dot(const Vector3& v) const {
        return m_vec.dot(v.m_vec);
    };

    // cross product computation method
    Vector3 Vector3::cross(const Vector3& v) const {
        return Vector3(m_vec.cross3(v.m_vec));
    };


    // vector rotation method
    // the rotation of the vector is centered at the origin of the global
    // reference frame
    // this method should not be used to rotate points with respect to a
    // specific frame
    // the axis is expected to be normalized
    // (see AffineTransformation class)
    void Vector3::rotate(const Real& angle, const Vector3& axis) {

        // normalized axis
        EigenVector3 v(axis[X], axis[Y], axis[Z]);
        v.normalize();

        // rotation transformation as quaternion
        // Eigen doc: for a single vector transformation quaternion is the
        // preferred representation
        const Eigen::Quaternion<Real> rot_q(Eigen::AngleAxis<Real>(angle, v));
        const EigenVector3 rv(rot_q*
                              EigenVector3(m_vec.x(),m_vec.y(),m_vec.z()));
        m_vec[0] = rv.x();
        m_vec[1] = rv.y();
        m_vec[2] = rv.z();

    };


    // array representation method
    Array<Real> Vector3::packed_values() const {
        return Array<Real>(3, m_vec.data());
    };

    // std::vector representation method
    std::vector<Real> Vector3::std_vector() const {
        return std::vector<Real>(m_vec.data(), m_vec.data()+3);
    };

    // move operator
    Vector3& Vector3::operator=(Vector3&& v) {
        // Eigen types are not movable at this point so we rely on the swap
        m_vec.swap(v.m_vec);
        return *this;
    };

    // add operator
    Vector3& Vector3::operator+=(const Vector3& v) {
        m_vec += v.m_vec;
        return *this;
    };

    // subtract operator
    Vector3& Vector3::operator-=(const Vector3& v) {
        m_vec -= v.m_vec;
        return *this;
    };

    // scalar multiplication operator
    Vector3& Vector3::operator*=(const Real& a) {
        m_vec *= a;
        return *this;
    };
    

    // addition operator
    Vector3 operator+(const Vector3& v1, const Vector3& v2) {
        return Vector3(v1.m_vec+v2.m_vec);
    };

    // subtraction operator
    Vector3 operator-(const Vector3& v1, const Vector3& v2) {
        return Vector3(v1.m_vec-v2.m_vec);
    };

    // scalar multiplication operator
    Vector3 operator*(const Vector3& v, const Real& a) {
        return Vector3(a*v.m_vec);
    };

    // scalar multiplication operator
    Vector3 operator*(const Real& a, const Vector3& v)  {
        return Vector3(a*v.m_vec);
    };


    // output operator
    // the output operator prints the vector as {X, Y, Z} without newline
	std::ostream& operator<<(std::ostream& os, const Vector3& v) {
		os << '{' << std::fixed << std::setprecision(REAL_WIDTH)
            << v.m_vec[X] << ", " << v.m_vec[Y] << ", " << v.m_vec[Z];
		os << '}';
		return os;
	};


    // private constructor from an Eigen vector object
    Vector3::Vector3(const EigenVector4& v) : m_vec(v) {};

    
}

