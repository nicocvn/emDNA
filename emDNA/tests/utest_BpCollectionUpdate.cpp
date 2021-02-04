// BpCollectionUpdate unit tests
// Nicolas Clauvelin


#include <emDNA.h>
#include "utest_BpCollectionUpdate.h"


namespace {


    // BpCollectionUpdateFromDofs
    TEST_F(BpCollectionUpdateTest, BpCollectionUpdateFromDofs) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(bp_set_1);

        VectorN newdofs(updated_dofs.size()*6, FLOAT_INIT);
        for (Size i=0; i<updated_dofs.size(); ++i)
            newdofs.set_slice(6*i, updated_dofs[i].inline_vector());

        FreeBpCollection coll;
        coll.set_bp_collection(bp_coll);
        coll.update_with_new_inline_free_dofs(newdofs);



        // update from dofs
//        bp_coll.update_collection_from_bp_step_dofs(updated_dofs);
        bp_coll = coll.bp_collection();

        // checking
        for (Size i=0; i<bp_coll.n_of_base_pairs(); ++i) {

            BasePair bp = bp_coll.base_pair(i);
            BasePair bp2 = bp_set_2[i];

            EXPECT_NEAR(FLOAT_INIT,
                        (bp.origin()-bp2.origin()).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (bp.axis(I)-bp2.axis(I)).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (bp.axis(J)-bp2.axis(J)).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (bp.axis(K)-bp2.axis(K)).norm(),
                        ZERO_TOL);

        };

    };


    // BpCollectionUpdateFromParams
    TEST_F(BpCollectionUpdateTest, BpCollectionUpdateFromParams) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(bp_set_1);

        // update from params
        bp_coll.update_collection_from_bp_step_params(updated_prms);

        // checking
        for (Size i=0; i<bp_coll.n_of_base_pairs(); ++i) {

            BasePair bp = bp_coll.base_pair(i);
            BasePair bp2 = bp_set_2[i];

            EXPECT_NEAR(FLOAT_INIT,
                        (bp.origin()-bp2.origin()).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (bp.axis(I)-bp2.axis(I)).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (bp.axis(J)-bp2.axis(J)).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (bp.axis(K)-bp2.axis(K)).norm(),
                        ZERO_TOL);
            
        };
        
    };
    

}
