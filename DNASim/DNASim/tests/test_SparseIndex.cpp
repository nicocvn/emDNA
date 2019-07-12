// test_SparseIndex class
// Nicolas Clauvelin


#include <SparseIndex.h>
#include <test_SparseIndex.h>


namespace {


    // GlobalTest
    TEST_F(SparseIndexTest, GlobalTest) {

        // index from 1 to 100
        SparseIndex index(1, 101);

        // flag
        index.flag_index(20-index.index_first_value());

        // flag range
        index.flag_index(30-index.index_first_value(),
                         40-index.index_first_value());

        // size checking
        EXPECT_EQ(Size(100), index.n_index_values());
        EXPECT_EQ(Size(11), index.n_flagged_indices());
        EXPECT_EQ(Size(89), index.n_clean_indices());

        // iterator checking
        SparseIndex::const_iterator it;
        for (it = index.begin_clean(); it != index.end_clean(); ++it) {
            EXPECT_TRUE(*it != 20);
            EXPECT_TRUE(*it < 30 || *it >= 40);
        };

    };


}

