// test_PointCloudKdTree class
// Nicolas Clauvelin


#ifndef DNASim_test_PointCloudKdTree_h
#define DNASim_test_PointCloudKdTree_h


#include <gtest/gtest.h>


namespace {


    // testing class for PointCloudKdTree
    class PointCloudKdTreeTest : public ::testing::Test {

	protected:

		PointCloudKdTreeTest() {};
		virtual ~PointCloudKdTreeTest() {};

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

#endif  // DNASim_test_PointCloudKdTree_h
