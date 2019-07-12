// test_MatrixN class
// Nicolas Clauvelin


#ifndef test_MatrixN_h
#define test_MatrixN_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for MatrixN
	class MatrixNTest : public ::testing::Test {
		
	protected:
		
		MatrixNTest() {}
		virtual ~MatrixNTest() {}
		
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


#endif	// test_MatrixN_h
