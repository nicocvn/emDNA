// test_PRN class
// Nicolas Clauvelin


#ifndef test_PRN_h
#define test_PRN_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for PRN
	class PRNTest : public ::testing::Test {
		
	protected:
		
		PRNTest() {}
		virtual ~PRNTest() {}
		
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


#endif	// test_PRN_h
