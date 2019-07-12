// SphereShape class
// Nicolas Clauvelin


// spheric shape implementation


#ifndef DNASim_SphereShape_h
#define DNASim_SphereShape_h


#include <Shape.h>


namespace DNASim {


    // pointer typedefs
    class SphereShape;


	class SphereShape final : public Shape {
		
		
	public:
		
		// constructors
		// default constructor creates a sphere of unit radius located at the
		// world origin
        // string constructor expects Sphere::position::radius
		SphereShape() = default;
		SphereShape(const SphereShape& sphere) = default;
        SphereShape(SphereShape&& sphere) = default;
        SphereShape(const Vector3& position, Real R);
        explicit SphereShape(const std::string& shape_data);
		~SphereShape() = default;
		
		// copy and move operators
		SphereShape& operator=(const SphereShape& sphere) = default;
        SphereShape& operator=(SphereShape&& sphere) = default;
		
		// cloning method
		Shape_Ptr clone() const override;
		
		// sphere radius accessor/modifier
		const Real&	radius() const;
		void set_radius(Real R);
		
		// inside point checking method
		bool is_point_inside(const Vector3& pt) const;
		
		// mathematica output method
		void mm_output(std::ostream& output,
					   const std::string& opts = "Blue") const;

		// ODE geom building method
#ifdef WITH_ODE_COLLISION
		dGeomID ode_geom() const override;
#endif	// WITH_ODE_COLLISION


	private:
		
		Real	m_radius;	// sphere radius
		
		
	};
	
	
}


#endif	// DNASim_SphereShape_h

