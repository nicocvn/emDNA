// Segment class
// Nicolas Clauvelin


// 3-dimensional segment implementation


#ifndef DNASim_Segment_h
#define DNASim_Segment_h


#include "maths/Vector3.h"


namespace DNASim {
	
	
	typedef Vector3 Vertex;


    class AffineTransformation;
	
	
	class Segment {	
		
		
	public:
		
		// constructors
		// default constructor creates a zero segment
		Segment() = default;
		Segment(const Vertex& v1, const Vertex& v2);
		Segment(const Segment& s) = default;
        Segment(Segment&& s) = default;
		~Segment() = default;
		
		// copy and move operators
		Segment& operator=(const Segment& s) = default;
        Segment& operator=(Segment&& s) = default;
		
		// point accessor/modifier
		inline const Vertex& first() const {
			return m_points.first;
		};
		inline const Vertex& second() const {
			return m_points.second;
		};
		inline Vertex& first() {
			return m_points.first;
		};
		inline Vertex& second() {
			return m_points.second;
		};
		
		// segment length accessor
		Real length() const;
		
		// segment tangent
		Vector3 tangent() const;
		
		// segment midpoint
		Vector3 midpoint() const;
		
		// geometry static methods
		// the turning angle is defined in [0;pi]
		static Real	turning_angle(const Segment& s1, const Segment& s2);
		static Vector3 binormal_vector(const Segment& s1, const Segment& s2);
		static Vector3 curvature_vector(const Segment& s1, const Segment& s2);
		static
        AffineTransformation parallel_transport_rotation(const Segment& s1,
                                                         const Segment& s2);
		
		// segment normal projection axis
		// compute the normal axis of the base segment pointing towards the edge
		// segment
		// based on the midpoints joining vector
		static Vector3 normal_projection_axis(const Segment& base,
                                              const Segment& edge);
		
		
	private:
		
		std::pair<Vertex,Vertex> m_points;
		
		
	};
	
	
}


#endif	// DNASim_Segment_h
