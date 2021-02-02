// test_Vector3 class
// Nicolas Clauvelin


#include "test_Vector3.h"
#include <DNASim.h>
using namespace DNASim;


namespace {
	
	
	// ConstructorFromString
	TEST_F(Vector3Test, ConstructorFromString) {

		Vector3 v(Real(2.457), Real(1.023), Real(-0.789));

		// test for components
		EXPECT_EQ(Real(2.457), v[X]);
		EXPECT_EQ(Real(1.023), v[Y]);
		EXPECT_EQ(Real(-0.789), v[Z]);
		
	};


    // MethodStdVector
    TEST_F(Vector3Test, MethodStdVector) {

		Vector3 v;
        v[X] = PRN_Real<RandomSeed::Yes>::generate();
        v[Y] = PRN_Real<RandomSeed::Yes>::generate();
        v[Z] = PRN_Real<RandomSeed::Yes>::generate();
        std::vector<Real> cv(v.std_vector());

        EXPECT_EQ(cv.size(), Size(3));
        EXPECT_EQ(v[X], cv[0]);
        EXPECT_EQ(v[Y], cv[1]);
        EXPECT_EQ(v[Z], cv[2]);
        
	};
	
	
	// MethodCross
	TEST_F(Vector3Test, MethodCross) {
		
		Vector3 u(Real(1.24), Real(2.35), Real(-0.58));
		Vector3 v(Real(-3.56), Real(1.85), Real(-1.25));
		Vector3 w(Real(-1.8645), Real(3.6148), Real(10.66));
		
		// test for components
        EXPECT_GE(ZERO_EPS, (w-u.cross(v)).norm());
		
	};
	
	
	// MethodRenormalizeAndIsUnitary
	TEST_F(Vector3Test, MethodNormalize) {
		
		Vector3 u(Real(1.24), Real(2.35), Real(-0.58));
		u.normalize();
		
		// test for components
		EXPECT_GE(ZERO_EPS, std::fabs(u.norm()-Real(1.0)));
		
	};


    // OperatorMultiplyWithMatrix3
	TEST_F(Vector3Test, OperatorMultiplyWithMatrix3) {

        Real values[9] = {Real(1.08), Real(-2.63), Real(-2.56), Real(0.17),
            Real(0.8), Real(0.67), Real(-0.96), Real(1.12), Real(-1.76)};
		Matrix3 m;
		m << Array<Real>(9, values);

        // input vector
		Vector3 x(Real(-1.34), Real(-1.03), Real(0.77));

        // result
		Vector3 res(-0.7095000000000005,
                    -0.5359,
                    -1.2224000000000002);

		Vector3 mx = m*x;

        EXPECT_GE(ZERO_EPS, (mx-res).norm());

	};


    // RotateMethod
    TEST_F(Vector3Test, RotateMethod) {

        // rotation data
		const Vector3 axis(Real(-0.55), Real(0.139), Real(.818));
        const Real angle = Real(2.74);

        // result
		const Vector3 res(-1.62087888860326829360662949506,
                          -0.273633880882030750462463645165,
                          -0.545780292529578590671593857709);

		Vector3 x(Real(1.0), Real(1.0), Real(1.0));
        x.rotate(angle, axis);
        EXPECT_GE(ZERO_EPS, (x-res).norm());
        
	};


    // EigensystemMethod
    TEST_F(Vector3Test, EigensystemMethod) {

        // test matrix
        Real values[9] = {Real(1.08), Real(-2.63), Real(-2.56), Real(0.17),
            Real(0.8), Real(0.67), Real(-0.96), Real(1.12), Real(-1.76)};
		Matrix3 m;
		m << Array<Real>(9, values);

        // eigenvalues
        Eigenvalues evalues = m.eigenvalues();

        // exact results
        Vector3 mm_res({
            1.8113947917907762, 0.844604461381112, -2.5359992531718882
        });

    };

	
}
