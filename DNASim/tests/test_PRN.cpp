// test_PRN class
// Nicolas Clauvelin


#include <PRN_Real.h>
#include <PRN_Integer.h>
#include <test_PRN.h>
using namespace DNASim;


namespace {

	
	// MethodGenerateWithBounds_Real
	TEST_F(PRNTest, MethodGenerateWithBounds_Real) {

		for (Size i=0; i<100; ++i) {
			Real f = PRN_Real<RandomSeed::Yes>::generate_with_bounds(Real(4.),
                                                                     Real(7.));
			bool chk = (Real(4.0)<=f) && (f<=Real(7.0));
			EXPECT_TRUE(chk);
		};
		
	};
	
	
	// MethodGenerateWithBounds_Integer
	TEST_F(PRNTest, MethodGenerateWithBounds_Integer) {
		
		for (Size i=0; i<100; ++i) {
			Real f = PRN_Integer<RandomSeed::Yes>::generate_with_bounds(-2, 1);
			bool chk = (-2<=f) && (f<=1);
			EXPECT_TRUE(chk);
		};
		
	};

	
}

