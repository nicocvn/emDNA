// BpCollectionGradient unit tests
// Nicolas Clauvelin


#include <emDNA.h>
#include "utest_BpCollectionGradient.h"


namespace {


    // FullGradientComputation
    TEST_F(BpCollectionGradientTest, FullGradientComputation) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_coll.set_collection_dummy_sequence();
        bp_coll.set_sequence_dependence_model("AnisoDNA_304");

        // index manager
        IndexManager idx_mgr;
        idx_mgr.set_n_of_bp_step(bp_coll.n_of_bp_steps());
        idx_mgr.set_frozen_steps_indexes(std::vector<SizePair>());

        // full gradient
        VectorN gradient =
        BpCollectionGradient::free_collection_gradient(bp_coll,
                                                       idx_mgr,
                                                       false);

        // checking
        for (Size i=0; i<mm_gradients.size(); ++i)
            for (Size j=0; j<6; ++j) {
                Size idx = idx_mgr.dof_global_coordinate(i, (BpStepDof)j);
                EXPECT_NEAR(FLOAT_INIT, gradient[idx]-mm_gradients[i][j],
                            ZERO_TOL);
            }

    };


    // ImposedLastBpGradientComputation
    TEST_F(BpCollectionGradientTest, ImposedLastBpGradientComputation) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_coll.set_collection_dummy_sequence();
        bp_coll.set_sequence_dependence_model("AnisoDNA_304");

        // index manager
        IndexManager idx_mgr;
        idx_mgr.set_n_of_bp_step(bp_coll.n_of_bp_steps());
        idx_mgr.set_frozen_steps_indexes(std::vector<SizePair>());

        // gradient
        VectorN gradient =
        BpCollectionGradient::EEDR_collection_gradient(bp_coll,
                                                       idx_mgr,
                                                       false);

        // checking
        for (Size i=0; i<mm_imp_gradients.size(); ++i)
            for (Size j=0; j<6; ++j) {
                Size idx = idx_mgr.dof_global_coordinate(i, (BpStepDof)j);
                EXPECT_NEAR(FLOAT_INIT, gradient[idx]-mm_imp_gradients[i][j],
                            ZERO_TOL);
            };

    };


}
