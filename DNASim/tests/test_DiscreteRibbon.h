// test_DiscreteRibon class
// Nicolas Clauvelin


#ifndef DNASim_test_DiscreteRibbon_h
#define DNASim_test_DiscreteRibbon_h


#include <gtest/gtest.h>
#include <DNASim_Includes.h>
using namespace DNASim;


namespace {
	
	
	// testing class for DiscreteRibbon
	class DiscreteRibbonTest : public ::testing::Test {
		
	protected:
		
		DiscreteRibbonTest() {}
		virtual ~DiscreteRibbonTest() {}
		
		// set up before test
		// called after the constructor
		virtual void SetUp() {};
		
		// clean up before test
		// called after the destructor
		virtual void TearDown() {};

		
	};
	
	
}


#endif	// DNASim_test_DiscreteRibbon_h
