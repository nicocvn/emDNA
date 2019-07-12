// test_Triad class
// Nicolas Clauvelin


#ifndef test_Triad_h
#define test_Triad_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for Triad
	class TriadTest : public ::testing::Test {
		
	protected:
		
		TriadTest() {}
		virtual ~TriadTest() {}
		
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


#endif	// test_Triad_h
