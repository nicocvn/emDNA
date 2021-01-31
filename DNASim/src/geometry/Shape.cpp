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


#ifdef WITH_ODE_COLLISION
	// ODE library init
	void Shape::ode_init() {
		// thread cleanup is manual
		dInitODE2(dInitFlagManualThreadCleanup);
	};
	
	// ODE library closing
	void Shape::ode_close() {
		dCloseODE();
	};
	
	// collision detection method
    bool Shape::ode_is_colliding(const Shape_Ptr__v& shapes_set_1,
                                 const Shape_Ptr__v& shapes_set_2) {

        // ODE init
		dAllocateODEDataForThread(0);

		// collision detection
#define MAX_CONTACTS 1

        // contact data array
		dContactGeom contact[MAX_CONTACTS];

        // loop over set 1
        for (auto it1 = shapes_set_1.begin(); it1 != shapes_set_1.end(); ++it1)
        {

            // first geom
			dGeomID g1 = (*it1)->ode_geom();

            // loop over set 2
            for (auto it2 = shapes_set_2.begin(); it2 != shapes_set_2.end();
                 ++it2) {

                // second geom
				dGeomID g2 = (*it2)->ode_geom();

                // collision check
				Integer coll = dCollide(g1,
										g2,
										MAX_CONTACTS,		// size of contact[]
										contact,
										sizeof(dContactGeom));

                // cleanup
				dGeomDestroy(g2);
				if (coll != 0) {
					dGeomDestroy(g1);
					dCleanupODEAllDataForThread();
					return true;
				};
			};

            // cleanup
			dGeomDestroy(g1);
            
        };
#undef MAX_CONTACTS

		// ODE finish
		dCleanupODEAllDataForThread();

		return false;


    };
#endif	// WITH_ODE_COLLISION

	
}
