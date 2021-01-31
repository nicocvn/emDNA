// Shape class
// Nicolas Clauvelin


#include "maths/Vector3.h"
#include "geometry/Shape.h"


namespace DNASim {


	// class constructor with position and frame initialization
	Shape::Shape(const Triad& f) : m_frame(f) {};


	// position accessor/modifier
	Vector3 Shape::position() const {
		return m_frame.origin();
	};
	void Shape::set_position(const Vector3& p) {
		m_frame.set_origin(p);
	};


	// shape frame accessor/modifier
	const Triad& Shape::frame() const {
		return m_frame;
	};
	void Shape::set_frame(const Triad& f) {
		m_frame = f;
	};

	
}
