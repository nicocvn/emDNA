// CylinderShape class
// Nicolas Clauvelin


// cylinder shape implementation


#ifndef DNASim_CylinderShape_h
#define DNASim_CylinderShape_h


#include "geometry/Shape.h"


namespace DNASim {


    // pointer typedefs
    class CylinderShape;

	
	class CylinderShape final : public Shape {
		
		
	public:
		
		// constructors
		// default constructor creates a cylinder of unit radius and height
		// located at the world origin and framed using the world coordinate
		// system
        // string constructor expects Cylinder::origin::X::Y::Z::radius::height
		CylinderShape() = default;
		CylinderShape(const CylinderShape& cylinder) = default;
        CylinderShape(CylinderShape&& cylinder) = default;
        CylinderShape(const Triad& frame, Real R, Real H);
        explicit CylinderShape(const std::string& shape_data);
		~CylinderShape() = default;
		
		// copy and move operators
		CylinderShape& operator=(const CylinderShape& cylinder) = default;
        CylinderShape& operator=(CylinderShape&& cylinder) = default;
		
		// cloning method
		Shape_Ptr clone() const override;
		
		// cylinder radius accessor/modifier
		const Real&	radius() const;
		void set_radius(Real R);
		
		// cylinder height accessor/modifier
		const Real&	height() const;
		void set_height(Real H);
		
		// inside point checking method
		bool is_point_inside(const Vector3& pt) const override;
		
		// mathematica output method
		void mm_output(std::ostream& output,
					   const std::string& opts = "Blue") const override;

		
	private:
		
		Real	m_radius;	// cylinder radius
		Real	m_height;	// cylinder height
		
		
	};
	
	
}


#endif	// DNASim_CylinderShape_h
