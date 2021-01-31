//// test_CollisionDetection class
//// Nicolas Clauvelin
//
//
//#include <Vector3.h>
//#include <Triad.h>
//#include <Shape.h>
//#include <SphereShape.h>
//#include <CylinderShape.h>
//#include <test_CollisionDetection.h>
//using namespace DNASim;
//
//
//namespace {
//
//
//    // SphereIsPointInsideMethod
//    TEST_F(CollisionDetectionTest, SphereIsPointInsideMethod) {
//
//        SphereShape s1(Vector3(), Real(5.0));
//
//        // points
//        Vector3 p1(FLOAT_INIT, 4.999999, FLOAT_INIT);
//        Vector3 p2(FLOAT_INIT, 5.0, FLOAT_INIT);
//        Vector3 p3(FLOAT_INIT, 5.000001, FLOAT_INIT);
//
//        // tests
//        EXPECT_TRUE(s1.is_point_inside(p1));
//        EXPECT_FALSE(s1.is_point_inside(p2));
//        EXPECT_FALSE(s1.is_point_inside(p3));
//
//    };
//
//
//    // CylinderIsPointInsideMethod
//    TEST_F(CollisionDetectionTest, CylinderIsPointInsideMethod) {
//
//        CylinderShape cyl(Triad(), 2.5, 4.0);
//        cyl.set_position(Vector3(3.0,-2.0,1.0));
//
//        // points
//        Vector3 p1(Vector3(3.0,-2.0,1.0)+Vector3(2.4999,FLOAT_INIT,FLOAT_INIT));
//        Vector3 p2(Vector3(3.0,-2.0,1.0)+Vector3(3.0,FLOAT_INIT,FLOAT_INIT));
//        Vector3 p3(Vector3(3.0,-2.0,1.0)+Vector3(3.0,7.0,FLOAT_INIT));
//
//        // tests
//        EXPECT_TRUE(cyl.is_point_inside(p1));
//        EXPECT_FALSE(cyl.is_point_inside(p2));
//        EXPECT_FALSE(cyl.is_point_inside(p3));
//
//    };
//
//
//	// SphereSphereCollision
//	TEST_F(CollisionDetectionTest, SphereSphereCollision) {
//
//		// init
//		Shape::ode_init();
//
//        // shapes
//        Shape_Ptr__v v1;
//        v1.push_back(Shape_Ptr(new SphereShape(Vector3(), Real(5.0))));
//        Shape_Ptr__v v2;
//        v2.push_back(Shape_Ptr(new SphereShape(Vector3(Real(-10.0),
//                                                       FLOAT_INIT,
//                                                       FLOAT_INIT),
//                                               Real(2.0))));
//
//        // no collision check
//        EXPECT_FALSE(Shape::ode_is_colliding(v1, v2));
//
//        // border case
//        v2[0]->set_position(Vector3(Real(-7.01), FLOAT_INIT, FLOAT_INIT));
//        EXPECT_FALSE(Shape::ode_is_colliding(v1, v2));
//		v2[0]->set_position(Vector3(Real(7.00), FLOAT_INIT, FLOAT_INIT));
//		EXPECT_TRUE(Shape::ode_is_colliding(v1, v2));
//
//		// close
//		Shape::ode_close();
//
//	};
//
//
//	// SphereCylinderCollision
//	TEST_F(CollisionDetectionTest, SphereCylinderCollision) {
//
//		// close
//		Shape::ode_init();
//
//        // shapes
//        Shape_Ptr__v v1;
//        v1.push_back(Shape_Ptr(new SphereShape(Vector3(Real(-5.0),
//                                                       FLOAT_INIT,
//                                                       FLOAT_INIT),
//                                               Real(1.0))));
//        Shape_Ptr__v v2;
//        v2.push_back(Shape_Ptr(new CylinderShape(Triad(),
//                                                 Real(1.0), Real(4.0))));
//
//		// no collision check
//        EXPECT_FALSE(Shape::ode_is_colliding(v1,v2));
//
//		// border case
//        v1[0]->set_position(Vector3(Real(-2.01), FLOAT_INIT, FLOAT_INIT));
//        EXPECT_FALSE(Shape::ode_is_colliding(v1,v2));
//		v1[0]->set_position(Vector3(Real(-1.99), FLOAT_INIT, FLOAT_INIT));
//		EXPECT_TRUE(Shape::ode_is_colliding(v1,v2));
//
//		// close
//		Shape::ode_close();
//
//	};
//
//
//	// CylinderCylinderCollision
//	TEST_F(CollisionDetectionTest, CylinderCylinderCollision) {
//
//		// close
//		Shape::ode_init();
//
//		// cylinder 1
//		Triad cyl1;
//		cyl1.set_origin(Vector3(Real(78.4208), Real(57.0192), Real(56.581)));
//		cyl1.set_axis(I, Vector3(Real(0.09686), Real(-0.83454), Real(0.54237)));
//		cyl1.set_axis(J, Vector3(Real(-0.95586), Real(-0.22987),
//                                 Real(-0.18299)));
//		cyl1.set_axis(K, Vector3(Real(0.27738), Real(-0.50071),
//                                 Real(-0.81997)));
//		Real cyl1R = Real(75.);
//		Real cyl1H = Real(43.);
//		CylinderShape cylinder1(cyl1, cyl1R, cyl1H);
//
//		// cylinder 2
//		Triad cyl2;
//		cyl2.set_origin(Vector3(Real(-16.8078), Real(57.9449),
//									   Real(39.9879)));
//		cyl2.set_axis(I, Vector3(FLOAT_INIT, FLOAT_INIT, Real(1)));
//		cyl2.set_axis(J, Vector3(Real(1), FLOAT_INIT, FLOAT_INIT));
//		cyl2.set_axis(K, Vector3(FLOAT_INIT, Real(1), FLOAT_INIT));
//		Real cyl2R = Real(55.);
//		Real cyl2H = Real(135.);
//		CylinderShape cylinder2(cyl2, cyl2R, cyl2H);
//
//		// geometry checking
//		Real odeR, odeH;
//		dGeomID cyl1_id = cylinder1.ode_geom();
//		dGeomCylinderGetParams(cyl1_id, &odeR, &odeH);
//		EXPECT_EQ(cyl1R, odeR);
//		EXPECT_EQ(cyl1H, odeH);
//		dGeomDestroy(cyl1_id);
//
//		// no collision check
//		Shape_Ptr__v v1;
//		Shape_Ptr__v v2;
//		v1.push_back(Shape_Ptr(new CylinderShape(cylinder1)));
//		v2.push_back(Shape_Ptr(new CylinderShape(cylinder2)));
//		EXPECT_TRUE(Shape::ode_is_colliding(v1, v2));
//
//		// close
//		Shape::ode_close();
//
//	};
//
//
//}
