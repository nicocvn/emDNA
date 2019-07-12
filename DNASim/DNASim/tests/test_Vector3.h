// test_Vector3 class
// Nicolas Clauvelin


#ifndef test_Vector3_h
#define test_Vector3_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for Vector3
	class Vector3Test : public ::testing::Test {
		
	protected:
		
		Vector3Test() {}
		virtual ~Vector3Test() {}
		
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


#endif	// test_Vector3_h
