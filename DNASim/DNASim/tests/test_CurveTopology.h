// test_CurveTopology class
// Nicolas Clauvelin


#ifndef test_CurveTopology_h
#define test_CurveTopology_h


#include <gtest/gtest.h>


namespace {
	
	
	// testing class for CurveTopology
	class CurveTopologyTest : public ::testing::Test {
		
	protected:
		
		CurveTopologyTest() {}
		virtual ~CurveTopologyTest() {}
		
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


#endif	// test_CurveTopology_h
