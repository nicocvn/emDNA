// test_MatrixN class
// Nicolas Clauvelin


#include <VectorN.h>
#include <MatrixN.h>
#include <test_MatrixN.h>
using namespace DNASim;


namespace {
	
	
	// MethodResize
	TEST_F(MatrixNTest, MethodResize) {
		
		MatrixN M(3,5);
		M(0,0) = Real(1);
		M(0,1) = Real(1);
		M.resize(2,3);
		
		EXPECT_EQ(M(0,0), Real(0));
		EXPECT_EQ((Size)3, M.n_cols());
		EXPECT_EQ((Size)2, M.n_rows());

	};
	
	
	// MethodTranspose
	TEST_F(MatrixNTest, MethodTranspose) {

        Real values[9] = {Real(-0.22), Real(1.4), Real(0.59), Real(1.8),
            Real(1.6), Real(0.56), Real(1.6), Real(1.9), Real(-0.28)};

		MatrixN M(3), N(3);
		M << Array<Real>(9, values);
		N = M;
		N.transpose();
		
		// check all entries
		for (Size i=0; i<3; ++i)
			for (Size j=0; j<3; ++j)
				EXPECT_GE(ZERO_TOL_SQ,(N(j,i)-M(i,j))*(N(j,i)-M(i,j)));
		
	};
	
	
	// SelfMatrixMultiplication
	TEST_F(MatrixNTest, SelfMatrixMultiplication) {

        Real values1[9] = {Real(-0.22), Real(1.4), Real(0.59), Real(1.8),
            Real(1.6), Real(0.56), Real(1.6), Real(1.9), Real(-0.28)};
        Real values2[9] = {Real(-0.66), Real(0.094), Real(-0.59), Real(-0.25),
			Real(-1.4), Real(0.94), Real(-1.8), Real(-1.9), Real(0.25)};

		MatrixN M(3);
		M << Array<Real>(9, values1);
		
		MatrixN N(3);
		N << Array<Real>(9, values2);
		
		M *= N;
		
		// pre-computed result
        Real results[9] = {Real(-1.2668), Real(-3.10168),
            Real(1.5932999999999997), Real(-2.596),
            Real(-3.1348), Real(0.582),
            Real(-1.0270000000000001),
            Real(-1.9775999999999998), Real(0.7719999999999998)};
		MatrixN R(3);
        R << Array<Real>(9, results);
		
		// check all entries
		for (Size i=0; i<3; ++i)
			for (Size j=0; j<3; ++j)
				EXPECT_GE(ZERO_TOL_SQ,(R(i,j)-M(i,j))*(R(i,j)-M(i,j)));
		
	};
	
	
	// SelfMatrixAddition
	TEST_F(MatrixNTest, SelfMatrixAddition) {
		
		Real values1[9] = {Real(-0.22), Real(1.4), Real(0.59), Real(1.8),
            Real(1.6), Real(0.56), Real(1.6), Real(1.9), Real(-0.28)};
        Real values2[9] = {Real(-0.66), Real(0.094), Real(-0.59), Real(-0.25),
			Real(-1.4), Real(0.94), Real(-1.8), Real(-1.9), Real(0.25)};

		MatrixN M(3);
		M << Array<Real>(9, values1);

		MatrixN N(3);
		N << Array<Real>(9, values2);
		
		M += N;

        // pre-computed result
        Real results[9] = {Real(-0.88), Real(1.494), Real(0.), Real(1.55),
            Real(0.2), Real(1.5), Real(-0.2), Real(0.), Real(-0.03)};
		MatrixN R(3);
        R << Array<Real>(9, results);

		// check all entries
		for (Size i=0; i<3; ++i)
			for (Size j=0; j<3; ++j)
				EXPECT_GE(ZERO_TOL_SQ,(R(i,j)-M(i,j))*(R(i,j)-M(i,j)));
		
	};
	
	
}

