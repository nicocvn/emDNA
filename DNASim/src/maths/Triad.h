// Triad class
// Nicolas Clauvelin


// direct orthonormal triad implementation
// a triad is defined as:
//  - an origin (3D point)
//  - an orthonormal set of three vectors

// see supporting notes for Euler ZYZ angles implementation details


#ifndef DNASim_Triad_h
#define DNASim_Triad_h


#include "maths/Matrix4.h"


// enum type for triad axis
enum TriadAxis {
    I=0,
    J=1,
    K=2
};


namespace DNASim {
	
	
	class Matrix3;
    class Vector3;
    class AffineTransformation;


	class Triad {
		

	public:
		
		// constructors
		// default constructor creates a direct unitary orthonormal frame
		// centered in (0,0,0)
		Triad();
        Triad(const Triad& f) = default;
        Triad(Triad&& f) = default;
		Triad(const Vector3& origin, const Matrix3& column_axis);
        Triad(const std::initializer_list<Vector3>& triad_list);
        explicit Triad(const Matrix4& m);
        explicit Triad(const std::string& s); // {origin, {axisI, axisJ, axisK}}
        ~Triad() = default;

        // copy and move operators
		Triad& operator=(const Triad& triad) = default;
        Triad& operator=(Triad&& triad) = default;

        // triad origin and axis accessors
		Vector3 origin() const;
		Vector3 axis(const TriadAxis& axis_index) const;
        Matrix3 column_axes() const;

        // triad origin and axis modifiers
		void set_origin(const Vector3& origin);
		void set_axis(const TriadAxis& axis_index, const Vector3& v);

        // triad matrix representation
		// returns a 4x4 augmented matrix containing axis and origin as columns
		// (the last row of the matrix is 0, 0, 0, 1)
		const Matrix4& matrix_representation() const;
        Matrix4& matrix_representation();

        // triad manipulation methods
		void normalize();
		void orthogonalize();
        void transform(const AffineTransformation& aff_t);

        // triad express methods
		// express_point_in return the coordinates in the current frame
		// express_point_out return the coordinates in the reference frame
		Vector3 express_point_in(const Vector3& v) const;
		Vector3 express_point_out(const Vector3& v) const;
        Vector3 express_vector_in(const Vector3& v) const;
		Vector3 express_vector_out(const Vector3& v) const;

		// transformation computation method
		// returns the transformation matrix f1->f2 as a 4x4 matrix
        // f2 = f1*transformation_matrix(f1,f2)
		static Matrix4 transformation_matrix(const Triad& f1,
                                             const Triad& f2);

		// Euler angles static methods
		// Euler angles are computed according to the ZYZ scheme
		// angles are defined for triad1->triad2
		static Vector3 euler_zyz_angles(const Triad& triad1,
                                        const Triad& triad2);
		static
        AffineTransformation euler_zyz_transformation(const Triad& triad1,
                                                      const Triad& triad2);
        static
        AffineTransformation euler_zyz_transformation(const Triad& triad,
                                                      const Vector3& angles);

		// output operator
		friend std::ostream& operator<<(std::ostream& os, const Triad& f);

		
	private:

        Matrix4 m_triad_mat;
		
		
	};
	
	
}


#endif	// DNASim_Triad_h
