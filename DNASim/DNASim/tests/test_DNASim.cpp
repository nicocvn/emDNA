// test_DNASim
// Nicolas Clauvelin


// entry point for DNASim tests based on googletest


#include <gtest/gtest.h>

#include <test_EnhancedString.h>
#include <test_PRN.h>
#include <test_Vector3.h>
#include <test_VectorN.h>
#include <test_MatrixN.h>
#include <test_AffineTransformation.h>
#include <test_Triad.h>
#include <test_CollisionDetection.h>
#include <test_StepParameters.h>
#include <test_StepArray.h>
#include <test_CurveTopology.h>
#include <test_DiscreteRibbon.h>


// main function
int main(int argc, char **argv) {

    try {

        // init
        ::testing::InitGoogleTest(&argc, argv);

        std::cout << std::fixed << std::setprecision(REAL_WIDTH);
        std::cerr << std::fixed << std::setprecision(REAL_WIDTH);

        // run all tests
        return RUN_ALL_TESTS();

    }
    catch (DNASim_ExitException& e) {
        exit(e._exit_code);
    };

};
