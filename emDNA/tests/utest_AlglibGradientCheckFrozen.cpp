// AlglibGradientCheckFrozen unit tests
// Nicolas Clauvelin


#include <FreeBpCollection.h>
#include <AnchoredBpCollection.h>
#include <ClampedBpCollection.h>
#include <AlglibBpCollectionMinimizer.h>
#include <Minimizer_AlglibLBFGS.h>
#include <utest_AlglibGradientCheckFrozen.h>


#define GRADCHECK_STEP_SIZE Real(1e-4)


namespace {


    // MinicircleClampedFrozen
    TEST_F(AlglibGradientCheckFrozen, MinicircleClampedFrozen) {

        // bp collection
        BpCollection bp_collection = create_minicircle(60);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("IdealDNA");

        // froze some steps
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(10,37));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        auto bp_collection_interface = std::make_shared<ClampedBpCollection>();
        bp_collection_interface->set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.1);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr = bp_collection_interface;
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface->
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibLBFGS minimizer_lbfgs;
        minimizer_lbfgs.set_minimizer_starting_point(free_dofs);
        minimizer_lbfgs.set_minimizer_settings(minim_settings);
        minimizer_lbfgs.set_minimizer_function(alglib_obj_function);
        Integer idx;
        AlglibMinResults minim_results =
        minimizer_lbfgs.
        perform_gradient_checking(GRADCHECK_STEP_SIZE,
                                  &idx,
                                  (void*)&callback_data);

        // checking
        // we allow for a single iteration so we should hit MAXIT
        // idx has to be -1 (no wrong component)
        EXPECT_EQ(MAXIT, minim_results._return_code);
        EXPECT_EQ(-1, idx);
        
    };


    // MinicircleClampedFrozenEnd
    TEST_F(AlglibGradientCheckFrozen, MinicircleClampedFrozenEnd) {

        // bp collection
        BpCollection bp_collection = create_minicircle(60);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("IdealDNA");

        // froze some steps
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(43,60));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        auto bp_collection_interface = std::make_shared<ClampedBpCollection>();
        bp_collection_interface->set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.1);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr = bp_collection_interface;
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface->
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibLBFGS minimizer_lbfgs;
        minimizer_lbfgs.set_minimizer_starting_point(free_dofs);
        minimizer_lbfgs.set_minimizer_settings(minim_settings);
        minimizer_lbfgs.set_minimizer_function(alglib_obj_function);
        Integer idx;
        AlglibMinResults minim_results =
        minimizer_lbfgs.
        perform_gradient_checking(GRADCHECK_STEP_SIZE,
                                  &idx,
                                  (void*)&callback_data);

        // checking
        // we allow for a single iteration so we should hit MAXIT
        // idx has to be -1 (no wrong component)
        EXPECT_EQ(MAXIT, minim_results._return_code);
        EXPECT_EQ(-1, idx);
        
    };


    // MinicircleClampedFrozenStart
    TEST_F(AlglibGradientCheckFrozen, MinicircleClampedFrozenStart) {

        // bp collection
        BpCollection bp_collection = create_minicircle(60);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("IdealDNA");

        // froze some steps
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(0,32));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        auto bp_collection_interface = std::make_shared<ClampedBpCollection>();
        bp_collection_interface->set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.1);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr = bp_collection_interface;
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface->
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibLBFGS minimizer_lbfgs;
        minimizer_lbfgs.set_minimizer_starting_point(free_dofs);
        minimizer_lbfgs.set_minimizer_settings(minim_settings);
        minimizer_lbfgs.set_minimizer_function(alglib_obj_function);
        Integer idx;
        AlglibMinResults minim_results =
        minimizer_lbfgs.
        perform_gradient_checking(GRADCHECK_STEP_SIZE,
                                  &idx,
                                  (void*)&callback_data);

        // checking
        // we allow for a single iteration so we should hit MAXIT
        // idx has to be -1 (no wrong component)
        EXPECT_EQ(MAXIT, minim_results._return_code);
        EXPECT_EQ(-1, idx);
        
    };


    // private method for minicircle guess
    BpCollection AlglibGradientCheckFrozen::create_minicircle(Size n) const {

        // step parameters
        std::vector<BpStepParams> params;
        params.reserve(n);

        // twist
        Real twist = Real(360)/n;
        for (Size i=0; i<n; ++i) {
            VectorN p(6, FLOAT_INIT);
            p[0] = twist*std::sin(34.3*DEG_2_RAD*i);
            p[1] = twist*std::cos(34.3*DEG_2_RAD*i);
            p[2] = 34.3;
            p[3] = FLOAT_INIT;
            p[4] = FLOAT_INIT;
            p[5] = Real(3.4);
            params.push_back(BpStepParams(p));
        };

        // create base pairs
        std::vector<BasePair> bps = BpStepParams::rebuild_bps(params,
                                                              BasePair());
        bps.pop_back();
        bps.push_back(bps.front());

        // bp collection
        return BpCollection::collection_from_base_pairs(bps);

    };


#undef GRADCHECK_STEP_SIZE


}

