// CylinderShape class
// Nicolas Clauvelin


#include <Vector3.h>
#include <EnhancedString.h>
#include <CylinderShape.h>


namespace DNASim {


	// class constructor with full initialization (position, frame, radius and
	// height)
	CylinderShape::CylinderShape(const Triad& cframe, Real R, Real H) :
	Shape(cframe), m_radius(R), m_height(H) {};
    
    
    // class constructor with initialization from a string
    CylinderShape::CylinderShape(const std::string& shape_data) :
    Shape(), m_radius(Real(1)), m_height(Real(1)) {
        
        // string tokenizing
        std::vector<std::string> data_str(EnhancedString::
                                          tokenize_string(shape_data, ':'));
        
        // check for proper shape
        DS_ASSERT(data_str[0] == "Cylinder",
                  "trying to construct a CylinderShape from a badly formatted"
                  " string: " + shape_data);
        
        // check for data completeness
        DS_ASSERT(data_str.size() == 7,
                  "trying to construct a CylinderShape from a badly formatted"
                  " string: " + shape_data);
        
        // cylinder triad
        Triad cylinder_triad;
        cylinder_triad.set_origin(Vector3(data_str[1]));
        cylinder_triad.set_axis(I, Vector3(data_str[2]));
        cylinder_triad.set_axis(J, Vector3(data_str[3]));
        cylinder_triad.set_axis(K, Vector3(data_str[4]));
        
        // radius
        Real cylinder_radius(EnhancedString::
                             convert_from_string<Real>(data_str[5]));
        
        // height
        Real cylinder_height(EnhancedString::
                             convert_from_string<Real>(data_str[6]));
        
        // shape settings
        set_frame(cylinder_triad);
        set_radius(cylinder_radius);
        set_height(cylinder_height);
        
    };


	// cloning method
	Shape_Ptr CylinderShape::clone() const {
		return Shape_Ptr(new CylinderShape(*this));
	};


	// cylinder radius accessor/modifier
	const Real& CylinderShape::radius() const {
		return m_radius;
	};
	void CylinderShape::set_radius(Real R) {
		m_radius = R;
	};

	// cylinder height accessor/modifier
	const Real& CylinderShape::height() const {
		return m_height;
	};
	void CylinderShape::set_height(Real H) {
		m_height = H;
	};


	// inside point checking method
	bool CylinderShape::is_point_inside(const Vector3& pt) const {
		
		// express point coordinates in the cylinder frame
		Vector3 localpt(frame().express_point_in(pt));
		
		if (localpt[X]*localpt[X]+localpt[Y]*localpt[Y] < m_radius*m_radius &&
			localpt[Z]*localpt[Z] < Real(0.25)*m_height*m_height)
			return true;

		return false;
		
	};


	// mathematica output method
	void CylinderShape::mm_output(std::ostream& output,
								  const std::string& opts) const {
        
		output << opts << ", Cylinder[{";
		output << position()-(m_height/Real(2))*frame().axis(K) << ",";
		output << position()+(m_height/Real(2))*frame().axis(K) << "}, ";
		output << m_radius << "]";
		
	};


#ifdef WITH_ODE_COLLISION
	// ODE geom building method
	dGeomID CylinderShape::ode_geom() const {
		
		// create cylinder geom
		dGeomID cylinder = dCreateCylinder(0, m_radius, m_height);
		
		// set position
		const Vector3& pos = position();
		dGeomSetPosition(cylinder, pos[X], pos[Y], pos[Z]);
		
		// set orientation
		const Matrix4& frame_mat = frame().matrix_representation();
		dMatrix3 orientation;
		orientation[0]	= frame_mat(0,0);
		orientation[1]	= frame_mat(0,1);
		orientation[2]	= frame_mat(0,2);
		orientation[4]	= frame_mat(1,0);
		orientation[5]	= frame_mat(1,1);
		orientation[6]	= frame_mat(1,2);
		orientation[8]	= frame_mat(2,0);
		orientation[9]	= frame_mat(2,1);
		orientation[10]	= frame_mat(2,2);
		dGeomSetRotation(cylinder, orientation);
		return cylinder;
		
	};
#endif	// WITH_ODE_COLLISION
	
	
}
