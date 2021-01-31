// test_DNASimArchive class
// Nicolas Clauvelin


#ifndef test_DNASimArchive_h
#define test_DNASimArchive_h


#include <gtest/gtest.h>
#include <DNASim.h>
using namespace DNASim;


namespace {


	// testing class for DNASimArchive
	class DNASimArchiveTest : public ::testing::Test {

	protected:

		DNASimArchiveTest() {};
		virtual ~DNASimArchiveTest() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {};

		// clean up before test
		// called after the destructor
		virtual void TearDown() {
            // clean archive file
            remove("test.cereal");
        };

		// declare objects that will be used by all tests
		// ...

        Vector3 v3 = Vector3(Real(2.457), Real(1.023), Real(-0.789));

        VectorN vN = VectorN({2.457, 1.023, -0.789, +17.23});

        Triad f = Triad();

        StepParameters p = StepParameters(VectorN({
            0.25, -0.17, 34.3, 0, 0, 3.4
        }));

	};


}


#endif	// test_DNASimArchive_h
