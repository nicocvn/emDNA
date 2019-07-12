// AlglibGradientCheckFrozen unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_AlglibGradientCheckFrozen_h
#define emDNA_utest_AlglibGradientCheckFrozen_h


#include <gtest/gtest.h>


namespace {


	// testing class for gradient computation on large system with multiple
    // frozen steps
	class AlglibGradientCheckFrozen : public ::testing::Test {

	protected:

		AlglibGradientCheckFrozen() {};
		virtual ~AlglibGradientCheckFrozen() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {};

        // clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests


        // private method for minicircle guess
        BpCollection create_minicircle(Size n) const;

	};
    
    
}


#endif  // emDNA_utest_AlglibGradientCheckFrozen_h
