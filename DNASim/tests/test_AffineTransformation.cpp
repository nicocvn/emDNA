// test_AffineTransformation class
// Nicolas Clauvelin



#include "test_AffineTransformation.h"
#include <DNASim.h>
using namespace DNASim;


namespace {
	
	
	// RotateVector3
	TEST_F(AffineTransformationTest, RotateVector3) {
		
		Vector3 v(Real(-3.56), Real(1.85), Real(-1.25));

		Vector3 axis(Real(1.01), Real(5.41), Real(-2.41));
        Real angle = Real(0.28);
        AffineTransformation t = AffineTransformation::rotation(angle,
                                                                axis);
		
		Vector3 res(Real(-3.517065453075459),
                    Real(2.2856824561942504),
                    Real(-0.2539809209939061));
        Vector3 rotv = t.transform_vector(v);

        EXPECT_GE(ZERO_TOL,(rotv-res).norm());
		
	};
	
	
	// RotateVector3MultipleRotations
	TEST_F(AffineTransformationTest, RotateVector3MultipleRotations) {
		
		Vector3 v(Real(-3.56), Real(1.85), Real(-1.25));

		Vector3 axis1(Real(1.01), Real(5.41), Real(-2.41));
		Vector3 axis2(Real(-0.25), Real(0.0), Real(-1.17));
        Real angle1 = Real(0.28);
        Real angle2 = Real(1.23);
        AffineTransformation t = AffineTransformation::rotation(angle1,
                                                                axis1);
        t.prepend(AffineTransformation::rotation(angle2, axis2));

		Vector3 res(Real(0.7943461573940103),
                    Real(3.95556136914705),
                    Real(-1.1752227181027664));
		Vector3 rotv = t.transform_vector(v);
		
        EXPECT_GE(ZERO_TOL,(rotv-res).norm());
		
	};


    // InverseTransformation
	TEST_F(AffineTransformationTest, InverseTransformation) {

        // original point
        Vector3 p0("{0.8053078315924664, 0.09940998343094476,"
                   "-0.8510009448829319}");

        Vector3 u("{-1, 2.3, 1.72}");
        Real angle = Real(.28);
        Vector3 axis(Real(1), Real(1), Real(1));
        axis.normalize();

        // transformation
        AffineTransformation t = AffineTransformation::translation(u);
        t.append(AffineTransformation::rotation(angle, axis));
        Vector3 p1 = t.transform_point(p0);

        // result
        Vector3 res("{-0.3769990752250466, 2.6605065357752347,"
                    "0.7902094095902906}");

        // inverse transformation
        AffineTransformation inv = t.inverse();
        Vector3 invp1 = inv.transform_point(p1);

        EXPECT_GE(ZERO_EPS,(p1-res).norm());
        EXPECT_GE(ZERO_EPS,(p0-invp1).norm());
		
	};
	
	
}
