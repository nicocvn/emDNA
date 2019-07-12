// Segment class
// Nicolas Clauvelin


#include <AffineTransformation.h>
#include <Segment.h>


namespace DNASim {
	
	
    // class constructor with vertices initialization
	Segment::Segment(const Vertex& v1, const Vertex& v2) : m_points(v1,v2) {};


	// segment length accessor
	Real Segment::length() const {
		return (second()-first()).norm();
	};
	
	
	// segment tangent
	Vector3 Segment::tangent() const {
        Vector3 v(second()-first());
        v.normalize();
		return v;
	};
	
	
	// segment midpoint
	Vector3 Segment::midpoint() const {
		return Real(0.5)*(second()+first());
	};
	
	
	// geometry
	Real Segment::turning_angle(const Segment& s1, const Segment& s2) {
		return std::acos(s1.tangent().dot(s2.tangent()));
	};
	Vector3 Segment::binormal_vector(const Segment& s1,
                                     const Segment& s2) {
        Vector3 kb = s1.tangent().cross(s2.tangent());
        kb.normalize();
		return kb;
	};
	Vector3 Segment::curvature_vector(const Segment& s1,
                                      const Segment& s2) {
		return s1.tangent().cross(s2.tangent());
	};
	AffineTransformation Segment::parallel_transport_rotation(const Segment& s1,
                                                              const Segment& s2)
    {
        Real angle = Segment::turning_angle(s1, s2);
        Vector3 axis(Segment::binormal_vector(s1, s2));
		return AffineTransformation::rotation(angle, axis);
	};
	
	
	// segment normal projection axis
	Vector3 Segment::normal_projection_axis(const Segment& base,
                                            const Segment& edge) {
		Vector3 j(edge.midpoint()-base.midpoint());
        Vector3 k(j-j.dot(base.tangent())*base.tangent());
        k.normalize();
		return k;
	};
	
	
}

