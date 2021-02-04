// AlglibMinimizer unit tests
// Nicolas Clauvelin


#include <emDNA.h>
#include "utest_AlglibMinimizer.h"


namespace {

    
    // ImposedLastBpGradientComputation
    TEST_F(AlglibMinimizerTest, ImposedLastBpGradientComputation) {

        // bp collection
        BpCollection bp_coll =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_coll.set_collection_dummy_sequence();
        bp_coll.set_sequence_dependence_model("IdealDNA_304");

        // bp collection interface
        ClampedBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_coll);

        // alglib conversion
        VectorN inline_dofs(bp_collection_interface.
                            bp_collection_inline_free_dofs());
        alglib::real_1d_array x(real_1d_array_with_size(inline_dofs.size()));
        VectorN_TO_real_1d_array(inline_dofs, x);

        // computation
        Real E(FLOAT_INIT);
        alglib::real_1d_array grad;
        grad.setlength(x.length());
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<ClampedBpCollection>(bp_collection_interface);
        callback_data._opts = opts;
        alglib_obj_function(x, E, grad, &callback_data);

        // partition
        std::vector<VectorN> gradients;
        for (Size i=0; i<bp_coll.n_of_bp_steps()-1; ++i) {
            VectorN g(6, FLOAT_INIT);
            g[0] = grad[6*(Integer)i];
            g[1] = grad[6*(Integer)i+1];
            g[2] = grad[6*(Integer)i+2];
            g[3] = grad[6*(Integer)i+3];
            g[4] = grad[6*(Integer)i+4];
            g[5] = grad[6*(Integer)i+5];
            gradients.push_back(g);
        };

        // checking
        for (Size i=0; i<gradients.size(); ++i)
            EXPECT_NEAR(FLOAT_INIT, (gradients[i]-mm_imp_gradients[i]).norm(),
                        ZERO_TOL);
        
    };


}
