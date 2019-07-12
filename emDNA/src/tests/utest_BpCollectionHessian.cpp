// BpCollectionHessian unit tests
// Nicolas Clauvelin


#include <IndexManager.h>
#include <BpCollection.h>
#include <HessianFunctions.h>
#include <BpCollectionHessian.h>
#include <utest_BpCollectionHessian.h>


namespace {


    // SingleStepHessian
    TEST_F(BpCollectionHessianTest, SingleStepHessianDiagonalTerm) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_bp_step_params(test_params, BasePair());
        bp_coll.set_collection_dummy_sequence();
        bp_coll.set_sequence_dependence_model("IdealDNA");

        // index manager
        IndexManager idx_mgr;
        idx_mgr.set_n_of_bp_step(bp_coll.n_of_bp_steps());
        idx_mgr.set_frozen_steps_indexes(std::vector<SizePair>());

        // single step Hessian
        const MatrixN diag_single_H =
        BpCollectionHessian::free_collection_Hessian(bp_coll, idx_mgr);
//        HessianFunctions::
//        single_step_Hessian_contrib(0, bp_coll, idx_mgr).full_matrix();

        // checking
        for (Size i=0; i<6; ++i) {
            for (Size j=0; j<6; ++j) {
                EXPECT_NEAR(mm_single_step_Hessian(i,j),
                            diag_single_H(i,j),
                            ZERO_TOL);
            };
        };

    };


    // FreeCollectionHessian
    TEST_F(BpCollectionHessianTest, FreeCollectionHessian) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_bp_step_params(test_params, BasePair());
        bp_coll.set_collection_sequence("ATCG");
        bp_coll.set_sequence_dependence_model("Olson1998");

        // index manager
        IndexManager idx_mgr;
        idx_mgr.set_n_of_bp_step(bp_coll.n_of_bp_steps());
        idx_mgr.set_frozen_steps_indexes(std::vector<SizePair>());

        // Hessian matrix
        const MatrixN H =
        BpCollectionHessian::free_collection_Hessian(bp_coll, idx_mgr);

        // checking
        for (Size i=0; i<3*6; ++i)
            for (Size j=0; j<3*6; ++j)
                EXPECT_NEAR(mm_free_Hessian(i,j),
                            H(i,j),
                            ZERO_TOL);

    };


}

