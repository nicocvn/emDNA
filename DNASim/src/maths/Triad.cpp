// Triad class
// Nicolas Clauvelin


#include "maths/Matrix3.h"
#include "maths/Vector3.h"
#include "maths/AffineTransformation.h"
#include "file_io/EnhancedString.h"
#include "maths/Triad.h"


namespace DNASim {


	// class default constructor
	// creates a triad locate in (0,0,0) with canonical axes
	Triad::Triad() : m_triad_mat(Matrix4::identity_matrix()) {
        orthogonalize();
    };


    // class constructor with origin and axis initialization
    Triad::Triad(const Vector3& origin, const Matrix3& column_axis) :
    m_triad_mat() {
        set_origin(origin);
        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                m_triad_mat(i,j) = column_axis(i,j);
    };


    // class constructor from a Vector3 initializer list
    Triad::Triad(const std::initializer_list<Vector3>& triad_list) {

        // size checking
        DS_ASSERT(triad_list.size()==4,
                  "wrong size for initializer list (4 is expected)");

        // setup
        set_origin(*(triad_list.begin()));
        set_axis(I, *(triad_list.begin()+1));
        set_axis(J, *(triad_list.begin()+2));
        set_axis(K, *(triad_list.begin()+3));

    };


    // class constructor with initialization from a 4x4 matrix
    Triad::Triad(const Matrix4& m) : m_triad_mat(m) {};


    // class constructor with initialization from a string
    Triad::Triad(const std::string& s) : m_triad_mat() {

        // string cleaning and tokenizing
        EnhancedString es(s);
        es.erase_blank_characters();
        es.erase_characters({'{','}'});
        std::vector<std::string> tokens;
        es.tokenize(tokens, ',');
        DS_ASSERT(tokens.size()==12,
                  "Triad(std::string) called with wrong size;"
				  " size=" << tokens.size());

        // origin and axes
        set_origin(Vector3(EnhancedString::convert_from_string<Real>(tokens[0]),
                           EnhancedString::convert_from_string<Real>(tokens[1]),
                           EnhancedString::convert_from_string<Real>(tokens[2]))
                   );
        set_axis(I,
                 Vector3(EnhancedString::convert_from_string<Real>(tokens[3]),
                         EnhancedString::convert_from_string<Real>(tokens[4]),
                         EnhancedString::convert_from_string<Real>(tokens[5]))
                 );
        set_axis(J,
                 Vector3(EnhancedString::convert_from_string<Real>(tokens[6]),
                         EnhancedString::convert_from_string<Real>(tokens[7]),
                         EnhancedString::convert_from_string<Real>(tokens[8]))
                 );
        set_axis(K,
                 Vector3(EnhancedString::convert_from_string<Real>(tokens[9]),
                         EnhancedString::convert_from_string<Real>(tokens[10]),
                         EnhancedString::convert_from_string<Real>(tokens[11]))
                 );

    };

    
    // triad origin and axis accessors
    Vector3 Triad::origin() const {
        return Vector3(m_triad_mat(0,3),
                       m_triad_mat(1,3),
                       m_triad_mat(2,3));
    };
    Vector3 Triad::axis(const TriadAxis& axis_index) const {
        return Vector3(m_triad_mat(0,axis_index),
                       m_triad_mat(1,axis_index),
                       m_triad_mat(2,axis_index));
    };
    Matrix3 Triad::column_axes() const {
        Real column_axes[9];
        column_axes[0] = m_triad_mat(0,0);
        column_axes[1] = m_triad_mat(0,1);
        column_axes[2] = m_triad_mat(0,2);
        column_axes[3] = m_triad_mat(1,0);
        column_axes[4] = m_triad_mat(1,1);
        column_axes[5] = m_triad_mat(1,2);
        column_axes[6] = m_triad_mat(2,0);
        column_axes[7] = m_triad_mat(2,1);
        column_axes[8] = m_triad_mat(2,2);
        Matrix3 m;
        m << Array<Real>(9, column_axes);
        return m;
    };


    // triad origin and axis modifiers
    void Triad::set_origin(const Vector3& origin) {
        m_triad_mat(0,3) = origin[X];
        m_triad_mat(1,3) = origin[Y];
        m_triad_mat(2,3) = origin[Z];
    };
    void Triad::set_axis(const TriadAxis& axis_index, const Vector3& v) {
        m_triad_mat(0,axis_index) = v[X];
        m_triad_mat(1,axis_index) = v[Y];
        m_triad_mat(2,axis_index) = v[Z];
    };


    // triad matrix representation
    // returns a 4x4 augmented matrix containing axis and origin as columns
    // (the last row of the matrix is 0, 0, 0, 1)
    const Matrix4& Triad::matrix_representation() const {
        return m_triad_mat;
    };
    Matrix4& Triad::matrix_representation() {
        return m_triad_mat;
    };


    // triad normalize method
	void Triad::normalize() {
        for (Size i=0; i<3; ++i) {
            Vector3 v = axis((TriadAxis)i);
            v.normalize();
            set_axis((TriadAxis)i, v);
        };
	};


	// triad orthogonalize method
	// based on Gram-Schmidt procedure
	void Triad::orthogonalize() {

        // current axis
        const Vector3 axis_I = axis(I);
        const Vector3 axis_J = axis(J);
        const Vector3 axis_K = axis(K);

		// build orthonormal axes
        Vector3 u1(axis_I);
        u1.normalize();             // new axis I
        Vector3 u2(axis_J-(u1.dot(axis_J))*u1);
        u2.normalize();             // new axis J
        Vector3 u3(axis_K-(u1.dot(axis_K))*u1-(u2.dot(axis_K))*u2);
        u3.normalize();             // new axis K

		// set triad axes
        m_triad_mat(0,I) = u1[X];
        m_triad_mat(1,I) = u1[Y];
        m_triad_mat(2,I) = u1[Z];
        m_triad_mat(0,J) = u2[X];
        m_triad_mat(1,J) = u2[Y];
        m_triad_mat(2,J) = u2[Z];
        m_triad_mat(0,K) = u3[X];
        m_triad_mat(1,K) = u3[Y];
        m_triad_mat(2,K) = u3[Z];
		
	};


    // frame transformation method
	void Triad::transform(const AffineTransformation& aff_t) {
        m_triad_mat = aff_t.matrix()*m_triad_mat;
	};


    // triad express methods
    // express_point_in return the coordinates in the current frame
    // express_point_out return the coordinates in the reference frame
    // in the following methods the _in methods use the inverse transform,
    // because we use row axis matrices
    Vector3 Triad::express_point_in(const Vector3& v) const {;
        const AffineTransformation affinv =
            AffineTransformation(m_triad_mat).inverse();
        return affinv.transform_point(v);
    };
    Vector3 Triad::express_point_out(const Vector3& v) const {
        const AffineTransformation aff_t(m_triad_mat);
        return aff_t.transform_point(v);
    };
    Vector3 Triad::express_vector_in(const Vector3& v) const {
        const AffineTransformation affinv =
            AffineTransformation(m_triad_mat).inverse();
        return affinv.transform_vector(v);
    };
    Vector3 Triad::express_vector_out(const Vector3& v) const {
        const AffineTransformation aff_t(m_triad_mat);
        return aff_t.transform_vector(v);
    };


	// transformation matrix computation method
	Matrix4 Triad::transformation_matrix(const Triad& f1, const Triad& f2) {

        // inverse affine transformation for f1
        const AffineTransformation inv_f1 =
            AffineTransformation(f1.matrix_representation()).inverse();

        // transformation matrix
        return inv_f1.matrix()*f2.matrix_representation();

	};


	// Euler angles ZYZ computation method
	// compute the euler angles corresponding to the rotation mapping triad1
	// onto triad2 within the ZYZ convention
	Vector3 Triad::euler_zyz_angles(const Triad& triad1,
                                    const Triad& triad2) {
		
		// needed scalar products
		Real kappacos = triad1.axis(K).dot(triad2.axis(K)); // D_33
		const Real alpha1 = triad1.axis(I).dot(triad2.axis(K));	// D_13
		const Real alpha2 = triad1.axis(J).dot(triad2.axis(K));	// D_23
		const Real lambda1	= triad1.axis(K).dot(triad2.axis(I)); // D_31
		const Real lambda2	= triad1.axis(K).dot(triad2.axis(J)); // D_32
		
		// check for 2D rotation case
		if (std::hypot(lambda1, lambda2) < ZERO_EPS) {

            const Real x = triad1.axis(I).dot(triad2.axis(I));
            const Real y = triad1.axis(J).dot(triad2.axis(I));
			
			// rotation angle
			const Real zeta = std::atan2(y, x);
			
			return Vector3(zeta/Real(2.0), Real(0.0), zeta/Real(2.0));
			
		};

        // safety net for kappacos
        if (kappacos > Real(1.0))
            kappacos = Real(1.0);
        if (kappacos < Real(-1.0))
            kappacos = Real(-1.0);

		// euler angles computation
		const Real kappa = std::acos(kappacos);
		const Real eta = std::atan2(lambda2, -lambda1);
		Real zeta = std::atan2(alpha2, alpha1);

		// smallest angle determination correction
		if (zeta+eta > F_PI)
			zeta -= Real(2)*F_PI;
		else if (zeta+eta < -F_PI)
			zeta += Real(2)*F_PI;

		return Vector3(zeta, kappa, eta);
		
	};
    AffineTransformation Triad::euler_zyz_transformation(const Triad& triad1,
                                                         const Triad& triad2) {

        // angles
        const Vector3 euler_angles = Triad::euler_zyz_angles(triad1, triad2);

        // axis
        const Vector3 Yaxis = triad1.axis(J);
        const Vector3 Zaxis = triad1.axis(K);

        // transformation
        AffineTransformation euler_t =
            AffineTransformation::rotation(euler_angles[X], Zaxis);
        euler_t.append(AffineTransformation::rotation(euler_angles[Y], Yaxis));
        euler_t.append(AffineTransformation::rotation(euler_angles[Z], Zaxis));

        return euler_t;

    };
    AffineTransformation Triad::euler_zyz_transformation(const Triad& triad,
                                                         const Vector3& angles)
    {

        // axis
        const Vector3 Yaxis = triad.axis(J);
        const Vector3 Zaxis = triad.axis(K);

        // transformation
        AffineTransformation euler_t =
            AffineTransformation::rotation(angles[X], Zaxis);
        euler_t.append(AffineTransformation::rotation(angles[Y], Yaxis));
        euler_t.append(AffineTransformation::rotation(angles[Z], Zaxis));

        return euler_t;

    };


    // friend output operator
	std::ostream& operator<<(std::ostream& os, const Triad& f) {
		os << '{';
		os << f.origin();
		os << ',';
		os << '{';
		os << f.axis(I) << ',' << f.axis(J) << ',' << f.axis(K);
		os << '}' << '}';
		return os;
	};


}
