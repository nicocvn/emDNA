// AlglibGradientCheck unit tests
// Nicolas Clauvelin


#include <emDNA.h>
#include "utest_AlglibGradientCheck.h"


#define GRADCHECK_STEP_SIZE Real(1e-8)


namespace {


    // MinCGGradientClampedCollection
    TEST_F(AlglibGradientCheckTest, MinCGGradientClampedCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        //bp_collection.set_sequence_dependence_model("TetramericModel"); // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // bp collection interface
        ClampedBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<ClampedBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibCG minimizer_cg;
        minimizer_cg.set_minimizer_starting_point(free_dofs);
        minimizer_cg.set_minimizer_settings(minim_settings);
        minimizer_cg.set_minimizer_function(alglib_obj_function);
        Integer idx;
        AlglibMinResults minim_results =
        minimizer_cg.perform_gradient_checking(GRADCHECK_STEP_SIZE,
                                               &idx,
                                               (void*)&callback_data);

        // checking
        // we allow for a single iteration so we should hit MAXIT
        // idx has to be -1 (no wrong component)
        EXPECT_EQ(MAXIT, minim_results._return_code);
        EXPECT_EQ(-1, idx);
        
    };


    // MinCGGradientAnchoredCollection
    TEST_F(AlglibGradientCheckTest, MinCGGradientAnchoredCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // bp collection interface
        AnchoredBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<AnchoredBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibCG minimizer_cg;
        minimizer_cg.set_minimizer_starting_point(free_dofs);
        minimizer_cg.set_minimizer_settings(minim_settings);
        minimizer_cg.set_minimizer_function(alglib_obj_function);
        Integer idx;
        AlglibMinResults minim_results =
        minimizer_cg.perform_gradient_checking(GRADCHECK_STEP_SIZE,
                                               &idx,
                                               (void*)&callback_data);

        // checking
        // we allow for a single iteration so we should hit MAXIT
        // idx has to be -1 (no wrong component)
        EXPECT_EQ(MAXIT, minim_results._return_code);
        EXPECT_EQ(-1, idx);
        
    };


    // MinLBFGSGradientClampedCollection
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientClampedCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // bp collection interface
        ClampedBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<ClampedBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
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


    // MinLBFGSGradientClampedFrozenCollection
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientClampedFrozenCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // froze middle step
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(3,3));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        ClampedBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<ClampedBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
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


    // MinLBFGSGradientClampedFrozenCollection
    TEST_F(AlglibGradientCheckTest,
           MinLBFGSGradientClampedFrozenCollectionElectrostatics) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // froze middle step
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(3,3));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        ClampedBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);
        bp_collection_interface.toggle_on_electrostatics();

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<ClampedBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
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


    // MinLBFGSGradientAnchoredCollection
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientAnchoredCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // bp collection interface
        AnchoredBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<AnchoredBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
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


    // MinLBFGSGradientAnchoredFrozenCollection
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientAnchoredFrozenCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // froze middle step
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(2,2));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        AnchoredBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<AnchoredBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
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


    // MinLBFGSGradientFreeCollection
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientFreeCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // bp collection interface
        FreeBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<FreeBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibLBFGS minimizer_lbfgs;
        minimizer_lbfgs.set_minimizer_starting_point(free_dofs);
        minimizer_lbfgs.set_minimizer_settings(minim_settings);
        minimizer_lbfgs.
        set_minimizer_function(alglib_obj_function);
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


    // MinLBFGSGradientFreeFrozenCheck
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientFreeFrozenCheck) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("AnisoDNA");

        // froze middle step
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(0,0));
        frozens.push_back(SizePair(2,2));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        FreeBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<FreeBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibLBFGS minimizer_lbfgs;
        minimizer_lbfgs.set_minimizer_starting_point(free_dofs);
        minimizer_lbfgs.set_minimizer_settings(minim_settings);
        minimizer_lbfgs.
        set_minimizer_function(alglib_obj_function);
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


    // MinLBFGSGradientPullCollection
    TEST_F(AlglibGradientCheckTest, MinLBFGSGradientPullCollection) {

        // bp collection
        BpCollection bp_collection =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bp_collection.set_collection_dummy_sequence();
        bp_collection.set_sequence_dependence_model("IdealDNA");

        // froze middle step
        std::vector<SizePair> frozens;
        frozens.push_back(SizePair(0,0));
        frozens.push_back(SizePair(2,2));
        bp_collection.set_frozen_steps_domains(frozens);

        // bp collection interface
        PullingBpCollection bp_collection_interface;
        bp_collection_interface.set_bp_collection(bp_collection);
        bp_collection_interface.set_pulling_force(Vector3(0,0,1));

        // minimization default settings
        AlglibMinSettings minim_settings;
        Real delta(1e-4);
        minim_settings._threshold_dx = Real(0.0);
        minim_settings._threshold_f = Real(0.0);
        minim_settings._threshold_g = delta;
        minim_settings._max_step_size = Real(0.0);
        minim_settings._max_iterations = 1;

        // callback data
        CallbackOptions opts;
        opts._callback_output = false;
        CallbackData callback_data;
        callback_data._bp_coll_ptr =
        std::make_shared<PullingBpCollection>(bp_collection_interface);
        callback_data._opts = opts;

        // minimization
        // do not use scalings when checking gradient
        VectorN free_dofs(bp_collection_interface.
                          bp_collection_inline_free_dofs());
        Minimizer_AlglibLBFGS minimizer_lbfgs;
        minimizer_lbfgs.set_minimizer_starting_point(free_dofs);
        minimizer_lbfgs.set_minimizer_settings(minim_settings);
        minimizer_lbfgs.
        set_minimizer_function(alglib_obj_function);
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


}


#undef GRADCHECK_STEP_SIZE
