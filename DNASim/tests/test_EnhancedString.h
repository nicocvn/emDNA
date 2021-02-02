// test_EnhancedString class
// Nicolas Clauvelin


#ifndef test_EnhancedString_h
#define test_EnhancedString_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for EnhancedString
	class EnhancedStringTest : public ::testing::Test {
		
	protected:
		
		EnhancedStringTest() {}
		virtual ~EnhancedStringTest() {}
		
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


#endif	// test_EnhancedString_h
