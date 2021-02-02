// CylinderShape class
// Nicolas Clauvelin


#include "maths/Vector3.h"
#include "file_io/EnhancedString.h"
#include "geometry/CylinderShape.h"


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

	
}
