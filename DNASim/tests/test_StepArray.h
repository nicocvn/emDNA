// test_StepArray class
// Nicolas Clauvelin


#ifndef test_StepArray_h
#define test_StepArray_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for StepArray
	class StepArrayTest : public ::testing::Test {
		
	protected:
		
		StepArrayTest() {}
		virtual ~StepArrayTest() {}
		
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


#endif	// test_StepArray_h
