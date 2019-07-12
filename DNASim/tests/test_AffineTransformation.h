// test_AffineTransformation class
// Nicolas Clauvelin


#ifndef test_AffineTransformation_h
#define test_AffineTransformation_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for AffineTransformation
	class AffineTransformationTest : public ::testing::Test {
		
	protected:
		
		AffineTransformationTest() {}
		virtual ~AffineTransformationTest() {}
		
		// set up before test
		// called after the constructor
		virtual void SetUp() {};
		
		// clean up before test
		// called after the destructor
		virtual void TearDown() {};
		
		// declare objects that will be used by all tests
		// ...
		
	};
	
	
}


#endif	// test_AffineTransformation_h
