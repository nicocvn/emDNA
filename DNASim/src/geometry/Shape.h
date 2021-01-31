// Shape class
// Nicolas Clauvelin


// abstract class for geometric shape implementation

// collision detection is based on ODE library


#ifndef DNASim_Shape_h
#define DNASim_Shape_h


#include "maths/Triad.h"


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
