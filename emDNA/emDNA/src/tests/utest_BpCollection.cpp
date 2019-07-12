// BpCollection unit tests
// Nicolas Clauvelin


#include <BpCollection.h>
#include <utest_BpCollection.h>


namespace  {


    // BpStepParamsComputation
    TEST_F(BpCollectionTest, BpStepParamsComputation) {

        // bp collection
        // step parameters are computed automatically
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(test_base_pairs);

        // check for step parameters values
        for (Size i=0; i<bp_coll.n_of_bp_steps(); ++i) {
            EXPECT_NEAR(test_bp_step_params[i].value(TILT),
                        bp_coll.bp_step_params(i).value(TILT),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(ROLL),
                        bp_coll.bp_step_params(i).value(ROLL),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(TWIST),
                        bp_coll.bp_step_params(i).value(TWIST),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(SHIFT),
                        bp_coll.bp_step_params(i).value(SHIFT),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(SLIDE),
                        bp_coll.bp_step_params(i).value(SLIDE),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(RISE),
                        bp_coll.bp_step_params(i).value(RISE),
                        ZERO_TOL);

        };

    };


    // BpStepDofsComputation
    TEST_F(BpCollectionTest, BpStepDofsComputation) {

        // bp collection
        // step parameters are computed automatically
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(test_base_pairs);

        // check for step parameters values
        for (Size i=0; i<bp_coll.n_of_bp_steps(); ++i) {
            EXPECT_NEAR(test_bp_step_dofs[i].value(TILTrad),
                        bp_coll.bp_step_dofs(i).value(TILTrad),
                        ZERO_EPS);
            EXPECT_NEAR(test_bp_step_dofs[i].value(ROLLrad),
                        bp_coll.bp_step_dofs(i).value(ROLLrad),
                        ZERO_EPS);
            EXPECT_NEAR(test_bp_step_dofs[i].value(TWISTrad),
                        bp_coll.bp_step_dofs(i).value(TWISTrad),
                        ZERO_EPS);
            EXPECT_NEAR(test_bp_step_dofs[i].value(R1),
                        bp_coll.bp_step_dofs(i).value(R1),
                        ZERO_EPS);
            EXPECT_NEAR(test_bp_step_dofs[i].value(R2),
                        bp_coll.bp_step_dofs(i).value(R2),
                        ZERO_EPS);
            EXPECT_NEAR(test_bp_step_dofs[i].value(R3),
                        bp_coll.bp_step_dofs(i).value(R3),
                        ZERO_EPS);
        };

        
    };


    // BpCollectionFromStepDofs
    TEST_F(BpCollectionTest, BpCollectionFromStepDofs) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_bp_step_dofs(test_bp_step_dofs,
                                                   BasePair());

        // check for step parameters values
        for (Size i=0; i<bp_coll.n_of_bp_steps(); ++i) {
            EXPECT_NEAR(test_bp_step_params[i].value(TILT),
                        bp_coll.bp_step_params(i).value(TILT),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(ROLL),
                        bp_coll.bp_step_params(i).value(ROLL),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(TWIST),
                        bp_coll.bp_step_params(i).value(TWIST),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(SHIFT),
                        bp_coll.bp_step_params(i).value(SHIFT),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(SLIDE),
                        bp_coll.bp_step_params(i).value(SLIDE),
                        ZERO_TOL);
            EXPECT_NEAR(test_bp_step_params[i].value(RISE),
                        bp_coll.bp_step_params(i).value(RISE),
                        ZERO_TOL);
        };

        // check for base pairs
        for (Size i=0; i<bp_coll.n_of_base_pairs(); ++i) {
            BasePair rbp = bp_coll.base_pair(i);
            EXPECT_NEAR(FLOAT_INIT,
                        (rbp.origin()-test_base_pairs[i].origin()).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (rbp.axis(I)-test_base_pairs[i].axis(I)).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (rbp.axis(J)-test_base_pairs[i].axis(J)).norm(),
                        ZERO_TOL);
            EXPECT_NEAR(FLOAT_INIT,
                        (rbp.axis(K)-test_base_pairs[i].axis(K)).norm(),
                        ZERO_TOL);
        };

    };


    // RebuildBpFromBpStepParams
    TEST_F(BpCollectionTest, RebuildBpFromBpStepParams) {

        BasePair rbp = test_bp_step_params[2].
            rebuild_step_last_base_pair(test_base_pairs[2]);

        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.origin()-test_base_pairs[3].origin()).norm(),
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.axis(I)-test_base_pairs[3].axis(I)).norm(),
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.axis(J)-test_base_pairs[3].axis(J)).norm(),
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.axis(K)-test_base_pairs[3].axis(K)).norm(),
                    ZERO_EPS);

    };


    // RebuildBpFromBpStepDofs
    TEST_F(BpCollectionTest, RebuildBpFromBpStepDofs) {

        BasePair rbp = test_bp_step_dofs[2].
            rebuild_step_last_base_pair(test_base_pairs[2]);

        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.origin()-test_base_pairs[3].origin()).norm(),
                    ZERO_TOL);
        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.axis(I)-test_base_pairs[3].axis(I)).norm(),
                    ZERO_TOL);
        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.axis(J)-test_base_pairs[3].axis(J)).norm(),
                    ZERO_TOL);
        EXPECT_NEAR(FLOAT_INIT,
                    (rbp.axis(K)-test_base_pairs[3].axis(K)).norm(),
                    ZERO_TOL);

    };


    // BpStepParamsAndBasePairDegenerated
    TEST_F(BpCollectionTest, BpStepParamsAndBasePairDegenerated) {

        BasePair bp1;
        BasePair bp2;

        // zero dofs
        BpStepParams p(bp1, bp2);
        EXPECT_NEAR(FLOAT_INIT, p.inline_vector().norm(), ZERO_TOL);

        // twist dofs
        bp2.transform(AffineTransformation::
                      rotation(Real(-0.27), Vector3(FLOAT_INIT, FLOAT_INIT, 1)));
        BpStepParams pbis(bp1, bp2);
        EXPECT_NEAR(Real(-0.27)*RAD_2_DEG, pbis.value(TWIST), ZERO_TOL);

        // rebuild bp
        BasePair bp = pbis.rebuild_step_last_base_pair(bp1);
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


    // BpStepDofsAndBasePairDegenerated
    TEST_F(BpCollectionTest, BpStepDofsAndBasePairDegenerated) {

        BasePair bp1;
        BasePair bp2;

        // zero dofs
        BpStepDofs dofs(bp1, bp2);
        EXPECT_NEAR(FLOAT_INIT, dofs.inline_vector().norm(), ZERO_EPS);

        // twist dofs
        bp2.transform(AffineTransformation::
                      rotation(Real(-0.27), Vector3(FLOAT_INIT, FLOAT_INIT, 1)));
        BpStepDofs dofsbis(bp1, bp2);
        EXPECT_NEAR(Real(-0.27), dofsbis.value(TWISTrad), ZERO_EPS);

        // rebuild bp
        BasePair bp = dofsbis.rebuild_step_last_base_pair(bp1);
        EXPECT_NEAR(FLOAT_INIT,
                    (bp.origin()-bp2.origin()).norm(),
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    (bp.axis(I)-bp2.axis(I)).norm(),
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    (bp.axis(J)-bp2.axis(J)).norm(),
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    (bp.axis(K)-bp2.axis(K)).norm(),
                    ZERO_EPS);

    };


    // BpCollectionListOfFreeIndices



    // BpCollectionRebuildWithFrozenSteps
    TEST_F(BpCollectionTest, BpCollectionRebuildWithFrozenSteps) {

            

    };


}
