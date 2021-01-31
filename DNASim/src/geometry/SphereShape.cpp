// SphereShape class
// Nicolas Clauvelin


#include "maths/Vector3.h"
#include "file_io/EnhancedString.h"
#include "geometry/SphereShape.h"


namespace DNASim {


	// class constructor with position and radius initialization
	SphereShape::SphereShape(const Vector3& pos, Real R) :
	Shape(Triad()), m_radius(R) {
		set_position(pos);
	};

    
    // class constructor with initialization from a string
    SphereShape::SphereShape(const std::string& shape_data) :
    Shape(Triad()), m_radius(Real(1)) {
        
        // string tokenizing
        std::vector<std::string> data_str(EnhancedString::
                                          tokenize_string(shape_data, ':'));
        
        // check for proper shape
        DS_ASSERT(data_str[0] == "Sphere",
                  "trying to construct a SphereShape from a badly formatted"
                  " string: " + shape_data);
        
        // check for data completeness
        DS_ASSERT(data_str.size() == 3,
                  "trying to construct a SphereShape from a badly formatted"
                  " string: " + shape_data);
        
        // set position
        set_position(Vector3(data_str[1]));
        
        // set radius
        set_radius(EnhancedString::convert_from_string<Real>(data_str[2]));
        
    };


	// cloning method
	Shape_Ptr SphereShape::clone() const {
		return Shape_Ptr(new SphereShape(*this));
	};


	// sphere radius accessor/modifier
	const Real& SphereShape::radius() const {
		return m_radius;
	};
	void SphereShape::set_radius(Real R) {
		m_radius = R;
	};


	// inside point checking method
	bool SphereShape::is_point_inside(const Vector3& pt) const {
		return ((pt-position()).norm() < m_radius) ? true : false;
	};


	// mathematica output method
	void SphereShape::mm_output(std::ostream& output,
								const std::string& opts) const {
		
		output << opts << ", Sphere[";
		output << position() << ", ";
		output << m_radius << "]";
		
	};


#ifdef WITH_ODE_COLLISION
	// ODE geom building method
	dGeomID SphereShape::ode_geom() const {
		
		// create sphere geom
		dGeomID sphere = dCreateSphere(0, m_radius);
		
		// set poisition
		const Vector3& pos = position();
		dGeomSetPosition(sphere, pos[X], pos[Y], pos[Z]);
		
		return sphere;
		
	};
#endif	// WITH_ODE_COLLISION
	
	
}

