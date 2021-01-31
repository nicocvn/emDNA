// emDNA_UnitTests
// Nicolas Clauvelin


#include <gtest/gtest.h>
#include <SequenceDepenceModels.h>


// entry point
GTEST_API_ int main(int argc, char **argv) {
    std::cout << ">> emDNA Unit Test <<\n\n";
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
};
