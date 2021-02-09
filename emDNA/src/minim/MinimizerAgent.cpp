// MinimizerAgent
// Nicolas Clauvelin


// bp collection interfaces headers
#include "DNA/FreeBpCollection.h"
#include "DNA/AnchoredBpCollection.h"
#include "DNA/ClampedBpCollection.h"
#include "DNA/PullingBpCollection.h"

// minimization headers
#include "minim/AlglibBpCollectionMinimizer.h"
#include "minim/MinimizerAgent.h"


MinimizationResults MinimizerAgent::
alglib_minimization(const std::shared_ptr<BpCollection_Interface>& bp_intf_ptr,
                    const CallbackOptions& callback_opts,
                    const AlgligMinSettings_Ptr& settings) {

    // minimization settigns
    AlglibMinSettings minim_settings(default_minim_settings());

    // override settings
    if (settings != nullptr)
        minim_settings = *settings;

    // minimization
    AlglibMinResults minim_results =
    AlglibBpCollectionMinimizer::
    optimize_bp_collection(bp_intf_ptr,
                           minim_settings,
                           LBFGS_ALGO,
                           callback_opts);

    // final results
    MinimizationResults final_res;
    final_res._optimized_bp_collection = bp_intf_ptr->bp_collection();
    final_res._n_iterations = minim_results._iterations;
    final_res._initial_E = minim_results._initial_function;
    final_res._final_E = minim_results._optimal_function;

    // free dofs gradient
    VectorN dofs_grad = minim_results._optimal_gradient;
    dofs_grad.diagonal_scale(bp_intf_ptr->bp_collection_free_dofs_scalings());
    final_res._free_dofs_gradient = dofs_grad;

    // full gradient
    final_res._elastic_gradient =
    bp_intf_ptr->bp_collection_elastic_energy_gradient();

    // energy contribs
    final_res._energy_contribs = bp_intf_ptr->energy_contribs();

    // return code string
    std::stringstream ss;
    ss << minim_results._return_code;
    final_res._return_code = ss.str();
    
    return final_res;

};


// default minim settings
AlglibMinSettings MinimizerAgent::default_minim_settings() {

    AlglibMinSettings minim_settings;

    Real delta = 1e-4;
    // 1 Jul 2019: Due to sizes of DNA used, alter max_iterations to work with sequences > 300 bp
    //minim_settings._max_iterations = 20000;
    minim_settings._max_iterations = 500000;
    minim_settings._threshold_dx = 0;
    minim_settings._threshold_f = 0;
    minim_settings._threshold_g = delta;
    minim_settings._max_step_size = Real(0.1);

    return minim_settings;
    
};
