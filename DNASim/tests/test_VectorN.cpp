// test_VectorN class
// Nicolas Clauvelin


#include "test_VectorN.h"
#include <DNASim.h>
using namespace DNASim;


namespace {
	
		
	// ConstructorFromString
	TEST_F(VectorNTest, ConstructorFromString) {
		
		const std::string s = "{2.457, 1.023, -0.789, +17.23}";
		VectorN v(s);;

		// test for size and components
		EXPECT_EQ(4, (Integer)v.size());
		EXPECT_EQ(Real(2.457), v[0]);
		EXPECT_EQ(Real(1.023), v[1]);
		EXPECT_EQ(Real(-0.789), v[2]);
		EXPECT_EQ(Real(17.23), v[3]);
		
	};


    // ConstructorFromStdVector
	TEST_F(VectorNTest, ConstructorFromStdVector) {

        std::vector<Real> src_v(5);
        for (Size i=0; i<5; ++i)
            src_v[i] = PRN_Real<RandomSeed::Yes>::generate();
        VectorN v(src_v);

		// test for size and components
		EXPECT_EQ(v.size(), Size(5));
        for (Size i=0; i<5; ++i)
            EXPECT_EQ(src_v[i], v[i]);

	};


    // MethodStdVector
    TEST_F(VectorNTest, MethodStdVector) {

		VectorN v(5, FLOAT_INIT);
        for (Size i=0; i<5; ++i)
            v[i] = PRN_Real<RandomSeed::Yes>::generate();
        std::vector<Real> cv(v.std_vector());

        EXPECT_EQ(cv.size(), Size(5));
        for (Size i=0; i<5; ++i)
            EXPECT_EQ(v[i], cv[i]);

	};
	
	
	// MethodNorm
	TEST_F(VectorNTest, MethodNorm) {
		
		VectorN v;
		EXPECT_EQ(0., v.norm());
		
	};
	
	
	// MethodAdd
	TEST_F(VectorNTest, MethodAdd) {
		
		VectorN v(4, 1.0);
		VectorN u(4, 1.25);
		
		v += Real(2.45)*u;
		
		EXPECT_EQ(Real(4.0625), v[0]);
		EXPECT_EQ(Real(4.0625), v[1]);
		EXPECT_EQ(Real(4.0625), v[2]);
		EXPECT_EQ(Real(4.0625), v[3]);
		
	};
	
	
	// MethodDiagonalScale
	TEST_F(VectorNTest, MethodDiagonalScale) {
		
		VectorN v(4, 1.0);
		VectorN u(4, 1.25);
		
		v.diagonal_scale(u);
		
		EXPECT_GE(ZERO_EPS, (v-u).norm());
		
	};
	
	
	// MethodDot
	TEST_F(VectorNTest, MethodDot) {
		
		VectorN v(4, 1.0);
		VectorN u(4, 1.25);
		
		Real d = v.dot(u);
		
		EXPECT_EQ(5.0, d);
		
	};


    // MethodSlice
	TEST_F(VectorNTest, MethodSlice) {

		VectorN v(4, 1.0);
        VectorN u(v.slice(2, 4));

		EXPECT_EQ(2, (Integer)u.size());
        EXPECT_EQ(1.0, u[0]);
        EXPECT_EQ(1.0, u[1]);

        u[0] = 2;
        u[1] = 3;
        v.set_slice(0, u);

        EXPECT_EQ(2, v[0]);
        EXPECT_EQ(3, v[1]);

	};
	

	// OperatorAdd
	TEST_F(VectorNTest, OperatorAdd) {
		
		VectorN u(3, 1.0);
		VectorN v(3, -2.0);
		
		VectorN w = v+u;
		
		// test for size and components
		EXPECT_EQ(3, (Integer)w.size());
		EXPECT_EQ(-1.0, w[0]);
		EXPECT_EQ(-1.0, w[1]);
		EXPECT_EQ(-1.0, w[2]);
		
	};
	
	
	// OperatorSubtract
	TEST_F(VectorNTest, OperatorSubtract) {
		
		VectorN u(3, 1.0);
		VectorN v(3, -2.0);
		
		VectorN w = v-u;
		
		// test for size and components
		EXPECT_EQ(3, (Integer)w.size());
		EXPECT_EQ(-3.0, w[0]);
		EXPECT_EQ(-3.0, w[1]);
		EXPECT_EQ(-3.0, w[2]);
		
	};
	
	
	// OperatorCopy
	TEST_F(VectorNTest, OperatorCopy) {
		
		VectorN u(3, 1.0);
		VectorN v(5, -2.0);
		
		u = v;
		
		// test for size and components
		EXPECT_GE(ZERO_EPS, (v-u).norm());
		
	};


    // OperatorMultiplyWithMatrixN
	TEST_F(VectorNTest, OperatorMultiplyWithMatrixN) {

        Real values[9] = {Real(1.08), Real(-2.63), Real(-2.56), Real(0.17),
            Real(0.8), Real(0.67), Real(-0.96), Real(1.12), Real(-1.76)};
		MatrixN m(3);
		m << Array<Real>(9, values);

        // vector input
		VectorN x(3, FLOAT_INIT);
        x[0] = -1.34;
        x[1] = -1.03;
        x[2] = 0.77;

        // result
        VectorN res(3, FLOAT_INIT);
        res[0] = -0.7095000000000005;
        res[1] = -0.5359;
        res[2] = -1.2224000000000002;

		VectorN mx = m*x;

		EXPECT_GE(ZERO_EPS, (mx-res).norm());

	};
	
}


