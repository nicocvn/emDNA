// Shape class
// Nicolas Clauvelin


// abstract class for geometric shape implementation

// collision detection is based on ODE library


#ifndef DNASim_Shape_h
#define DNASim_Shape_h


#include <ODE_Includes.h>
#include <Triad.h>


namespace DNASim {


	class Vector3;

    // type renaming
    class Shape;
    using Shape_Ptr = std::unique_ptr<Shape>;
    using Shape_Ptr__v = std::vector<Shape_Ptr>;
    

    class Shape {

		
	public:
		
		// virtual destructor
        virtual ~Shape() = default;
		
		// cloning method
        virtual Shape_Ptr clone() const =0;
		
		// position accessor/modifier
		Vector3 position() const;
		void set_position(const Vector3& p);
		
		// shape frame accessor/modifier
		const Triad& frame() const;
		void set_frame(const Triad& f);
		
		// inside point checking method
		virtual bool is_point_inside(const Vector3& pt) const =0;
		
		// mathematica output method
		virtual void mm_output(std::ostream& output,
							   const std::string& opts = "Blue") const =0;

#ifdef WITH_ODE_COLLISION
		// if using ODE collision detection:
		// - Shape::ode_init should be called before any collision detection
		// - Shape::ode_close should be called when no more collision detection
		// will be performed
		// - should be thread-safe
		
		// ODE library init
		static void ode_init();
		
		// ODE library closing
		static void ode_close();
		
		// ODE geom building method
		virtual dGeomID ode_geom() const =0;
		
		// collision detection method
		// return true if collision is detected
        static bool ode_is_colliding(const Shape_Ptr__v& shapes_set_1,
                                     const Shape_Ptr__v& shapes_set_2);
#endif	// WITH_ODE_COLLISION
		
		
	protected:
		
		// constructors
		// default constructor creates a shape located at the world origin and
		// framed using the world coordinate system
		Shape() = default;
		Shape(const Shape& shape) = default;
        Shape(Shape&& shape) = default;
        explicit Shape(const Triad& frame);
		
		// copy and move operators
		Shape& operator=(const Shape& shape) = default;;
        Shape& operator=(Shape&& shape) = default;;
		
		
	private:
		
		Triad m_frame; // shape frame
		
		
	};
	
	
}


#endif	// DNASim_Shape_h
