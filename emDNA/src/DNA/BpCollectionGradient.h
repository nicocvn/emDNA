// BpCollectionGradient functions
// Nicolas Clauvelin


// implementation of gradient computations


// notes about gradient computation:
// - every gradients are computed as follow (the order matters):
//      1. compute the gradient with respect to all dofs
//      2. add contributions for the boundary conditions
//      3. add contributions for the frozen steps
//      4. filter the resulting gradient to select the free steps only
// - the computations of the different contributions are implemented as pencil
//   functions


#ifndef emDNA_BpCollectionGradient_h
#define emDNA_BpCollectionGradient_h


#include "DNA/StepGradient.h"
class BpStepParams;
class BpCollection;
class IndexManager;


namespace BpCollectionGradient {


    // gradient computation functions
    VectorN free_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr,
                                     bool elec_flag);
    VectorN EEDR_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr,
                                     bool elec_flag);
    VectorN EED_collection_gradient(const BpCollection& bp_collection,
                                    const IndexManager& idx_mgr,
                                    bool elec_flag);
    VectorN pulling_collection_gradient(const BpCollection& bp_collection,
                                        const IndexManager& idx_mgr,
                                        const Vector3& pulling_force,
                                        bool elec_flag);

    // single step parameters gradient
    StepGradient single_step_parameters_gradient(const BpStepParams& p,
                                                 const BpStepParams& p0,
                                                 const MatrixN& elastic_moduli);


    // following code is in an unnamed namespace
    namespace {


    // collection step parameters gradient
    // this gradient is computed over all free steps
    // (the step parameters gradient of a frozen step is zero by definition)
    StepGradientVec step_parameters_gradient(const BpCollection& bp_collection,
                                             const IndexManager& idx_mgr);

    // collection dofs gradient
    // this gradient is computed over all steps
    // this gradient is used as a basis for the computation of other gradients
    StepGradientVec dofs_gradient(const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr);

    // pencil function for EEDR contributions to the gradients
    // those contributions are computed over all steps
    void add_eedr_contributions(StepGradientVec& step_grads,
                                const BpCollection& bp_collection,
                                const IndexManager& idx_mgr);

    // pencil function for EED contributions to the gradients
    // those contributions are computed over all steps
    void add_eed_contributions(StepGradientVec& step_grads,
                               const BpCollection& bp_collection,
                               const IndexManager& idx_mgr);

    // pencil function for frozen contributions to the gradients
    // those contributions are computed over all steps
    void add_frozen_contributions(StepGradientVec& step_grads,
                                  const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr);

    // pencil function for electrostatics contributions to the gradients
    // those contributions are computed over all steps
    void add_electrostatics_contributions(StepGradientVec& step_grads,
                                          const BpCollection& bp_collection);


    }


}


#endif  // emDNA_BpCollectionGradient_h
