// test_VectorN class
// Nicolas Clauvelin


#ifndef test_VectorN_h
#define test_VectorN_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for VectorN
	class VectorNTest : public ::testing::Test {
		
	protected:
		
		VectorNTest() {}
		virtual ~VectorNTest() {}
		
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


#endif	// test_Vector_h

