// BpCollectionGradient functions
// Nicolas Clauvelin


#include <DNAElectrostaticsParams.h>
#include <IndexManager.h>
#include <BpCollection.h>
#include <BpGeometryFunctions.h>
#include <BpCollectionGeometry.h>
#include <BpCollectionGradient.h>


namespace BpCollectionGradient {


    // gradient computation function for free collection
    VectorN free_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr,
                                     bool elec_flag) {

        // dofs gradient
        StepGradientVec dofs_grad = dofs_gradient(bp_collection, idx_mgr);

        // electrostatics interactions
        if (elec_flag)
            add_electrostatics_contributions(dofs_grad, bp_collection);

        // frozen contributions
        if (idx_mgr.n_of_frozen_steps() != 0)
            add_frozen_contributions(dofs_grad, bp_collection, idx_mgr);

        // filtering
        const StepGradientVec final_grads =
        filter_free_dofs_grads(dofs_grad, idx_mgr);

        // trimming
        return flatten_StepGradientVec(final_grads);

    };


    // gradient computation function for collection with EEDR boundary
    // conditions
    VectorN EEDR_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr,
                                     bool elec_flag) {

        // dofs gradient
        StepGradientVec dofs_grad = dofs_gradient(bp_collection, idx_mgr);

        // electrostatics interactions
        if (elec_flag)
            add_electrostatics_contributions(dofs_grad, bp_collection);

        // EEDR contributions
        add_eedr_contributions(dofs_grad, bp_collection, idx_mgr);

        // frozen contributions
        if (idx_mgr.n_of_frozen_steps() != 0)
            add_frozen_contributions(dofs_grad, bp_collection, idx_mgr);

        // filtering
        const StepGradientVec final_grads =
        filter_free_dofs_grads(dofs_grad, idx_mgr);

        // trimming
        VectorN grad_vec = flatten_StepGradientVec(final_grads);
        return grad_vec.slice(0,
                              idx_mgr.n_of_free_dofs_variables()-
                              BoundaryConditionsDofsTrimming::EEDRCollection);

    };


    // gradient computation function for collection with EED boundary
    // conditions
    VectorN EED_collection_gradient(const BpCollection& bp_collection,
                                    const IndexManager& idx_mgr,
                                    bool elec_flag) {

        // dofs gradient
        StepGradientVec dofs_grad = dofs_gradient(bp_collection, idx_mgr);

        // electrostatics interactions
        if (elec_flag)
            add_electrostatics_contributions(dofs_grad, bp_collection);

        // EED contributions
        add_eed_contributions(dofs_grad, bp_collection, idx_mgr);

        // frozen contributions
        if (idx_mgr.n_of_frozen_steps() != 0)
            add_frozen_contributions(dofs_grad, bp_collection, idx_mgr);

        // filtering
        const StepGradientVec final_grads =
        filter_free_dofs_grads(dofs_grad, idx_mgr);

        // trimming
        VectorN grad_vec = flatten_StepGradientVec(final_grads);
        return grad_vec.slice(0,
                              idx_mgr.n_of_free_dofs_variables()-
                              BoundaryConditionsDofsTrimming::EEDCollection);

    };


    // gradient computation function for collection with pulling force on the
    // last bp
    VectorN pulling_collection_gradient(const BpCollection& bp_collection,
                                        const IndexManager& idx_mgr,
                                        const Vector3& pulling_force,
                                        bool elec_flag) {

        // dofs gradient
        StepGradientVec dofs_grad = dofs_gradient(bp_collection, idx_mgr);

        // electrostatics interactions
        if (elec_flag)
            add_electrostatics_contributions(dofs_grad, bp_collection);

        // add contributions from the pulling force
        for (Size i=0; i<bp_collection.n_of_bp_steps(); ++i) {
            dofs_grad[i]._translation[X] -=
            pulling_force[X]*emDNAConstants::PullingForceScaling;
            dofs_grad[i]._translation[Y] -=
            pulling_force[Y]*emDNAConstants::PullingForceScaling;
            dofs_grad[i]._translation[Z] -=
            pulling_force[Z]*emDNAConstants::PullingForceScaling;
        };

        // frozen contributions
        if (idx_mgr.n_of_frozen_steps() != 0)
            add_frozen_contributions(dofs_grad, bp_collection, idx_mgr);

        // filtering
        const StepGradientVec final_grads =
        filter_free_dofs_grads(dofs_grad, idx_mgr);

        return flatten_StepGradientVec(final_grads);

    };


    // single step parameters gradient
    StepGradient single_step_parameters_gradient(const BpStepParams& p,
                                                 const BpStepParams& p0,
                                                 const MatrixN& elastic_moduli)
    {

        // symmetrix elastic moduli matrix
        MatrixN symmetrize_fmat = elastic_moduli;
        symmetrize_fmat.transpose();
        symmetrize_fmat += elastic_moduli;

        // gradient
        VectorN g(Real(0.5)*(symmetrize_fmat*(p.inline_vector()-
                                              p0.inline_vector())));

        // gradient structure
        return StepGradient({
            { g[0]*RAD_2_DEG, g[1]*RAD_2_DEG, g[2]*RAD_2_DEG },
            { g[3], g[4], g[5] }
        });
        
    };


// following code is in an unnamed namespace
namespace {


    // bp collection step parameters gradient
    // this vector is computed over all free steps
    StepGradientVec step_parameters_gradient(const BpCollection& bp_collection,
                                             const IndexManager& idx_mgr) {

        // free steps iterators
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // container
        std::vector<StepGradient> grads;
        grads.reserve(idx_mgr.n_of_free_steps());

        // loop over free steps
        for (auto free_it = free_begin; free_it != free_end; ++free_it) {

            // step collection index
            const Size& coll_idx = free_it->collection_index();

            // step data
            const BpStepParams& p =
            bp_collection.bp_step_params(coll_idx);
            const BpStepParams& p0 =
            bp_collection.bp_step_intrinsic_parameters(coll_idx);
            const MatrixN& fmat =
            bp_collection.bp_step_force_constants(coll_idx);

            // single step gradient
            grads.push_back(std::
                            move(single_step_parameters_gradient(p, p0, fmat)));

        };

        return grads;

    };


    // collection dofs gradient
    // this gradient is computed over all steps
    StepGradientVec dofs_gradient(const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr) {

        // size
        const Size n_steps = bp_collection.n_of_bp_steps();

        // step parameters gradient
        // this gradient only contains free steps contributions
        const StepGradientVec sp_grad =
        std::move(step_parameters_gradient(bp_collection,
                                           idx_mgr));

        // free steps iterators
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // dofs gradient container
        StepGradientVec dofs_grad(n_steps, StepGradient());

        // single step contributions
        // we loop over all free steps
        for (auto free_it = free_begin; free_it != free_end; ++free_it) {

            // collection index
            const Size& coll_idx = free_it->collection_index();

            // free step index
            const Size& free_idx = free_it->local_index();

            // matrices R and T
            const Matrix3 R_mat =
            BpCollectionGeometry::jacobian_matrix_R(coll_idx, bp_collection);
            const Matrix3 T_mat =
            BpCollectionGeometry::jacobian_matrix_T(coll_idx, bp_collection);

            // local terms
            // this corresponds to the single step dofs gradient
            StepGradient g = {
                sp_grad[free_idx]._rotation +
                    R_mat*sp_grad[free_idx]._translation,
                T_mat*sp_grad[free_idx]._translation
            };

            dofs_grad[coll_idx] = std::move(g);

        };

        // coupling contributions
        // we loop over all steps and check for free downward steps
        for (Size k=0; k<n_steps; ++k) {

            // we iterate over all free steps
            for (auto free_it = free_begin; free_it != free_end; ++free_it) {

                // collection index
                const Size& coll_idx = free_it->collection_index();

                // free step index
                const Size& free_idx = free_it->local_index();

                // check if this is a downward step
                // if yes, we compute the coupling contribution
                if (coll_idx > k) {

                    // matrix U
                    const Matrix3 U_mat =
                    BpCollectionGeometry::jacobian_matrix_U(k, coll_idx,
                                                            bp_collection);

                    // add the contribution to the gradient
                    dofs_grad[k]._rotation +=
                    U_mat*sp_grad[free_idx]._translation;

                };

            };

        };

        return dofs_grad;

    };


    // pencil function for EEDR contributions to the gradients
    // those contributions are computed over all steps
    void add_eedr_contributions(StepGradientVec& step_grads,
                                const BpCollection& bp_collection,
                                const IndexManager& idx_mgr) {

        // size
        const Size n_steps = bp_collection.n_of_bp_steps();

        // size checking
        DS_ASSERT(step_grads.size()==n_steps,
                  "wrong size for the step gradients;"
                  "current size is " +
                  EnhancedString::convert_to_string(step_grads.size()));

        // step gradients copy
        const StepGradientVec input_grads = step_grads;

        // bc step index
        const Size& bc_index = idx_mgr.bc_reduced_step().collection_index();

        // loop over all steps preceeding the bc step
        for (Size i=0; i<bc_index; ++i) {

            // reduction matrix
            Matrix3 Kmat =
            BpCollectionGeometry::jacobian_matrix_EER(i, bc_index,
                                                      bp_collection);
            Kmat.transpose();

            // bc contribution
            step_grads[i]._rotation +=
            Real(-1.0)*Kmat*input_grads[bc_index]._rotation;
            step_grads[i]._translation +=
            Real(-1.0)*input_grads[bc_index]._translation;
            
        };

    };

    // pencil function for EED contributions to the gradients
    // those contributions are computed over all steps
    void add_eed_contributions(StepGradientVec& step_grads,
                               const BpCollection& bp_collection,
                               const IndexManager& idx_mgr) {

        // size
        const Size n_steps = bp_collection.n_of_bp_steps();

        // size checking
        DS_ASSERT(step_grads.size()==n_steps,
                  "wrong size for the step gradients;"
                  "current size is " +
                  EnhancedString::convert_to_string(step_grads.size()));

        // step gradients copy
        const StepGradientVec input_grads = step_grads;

        // loop over all steps preceeding the bc step
        const Size& bc_index =
        idx_mgr.bc_reduced_step().collection_index();
        for (Size i=0; i<bc_index; ++i) {

            // bc contribution
            step_grads[i]._translation +=
            Real(-1.0)*input_grads[bc_index]._translation;
            
        };

    };


    // pencil function for frozen contributions to the gradients
    // those contributions are computed over all steps
    void add_frozen_contributions(StepGradientVec& step_grads,
                                  const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr) {

        // size
        const Size n_steps = bp_collection.n_of_bp_steps();

        // size checking
        DS_ASSERT(step_grads.size()==n_steps,
                  "wrong size for the step gradients;"
                  "current size is " +
                  EnhancedString::convert_to_string(step_grads.size()));

        // frozen map computation

        // step gradients copy
        const StepGradientVec input_grads = step_grads;

        // loop over all frozen steps
        const auto frz_begin = idx_mgr.frozen_step_indexes_begin();
        const auto frz_end = idx_mgr.frozen_step_indexes_end();
        for (auto frz_it = frz_begin; frz_it != frz_end; ++frz_it) {

            // collection index
            const Size& coll_idx = frz_it->collection_index();

            // frozen step dofs gradient
            const StepGradient& frz_grad = input_grads[coll_idx];

            // loop on preceeding steps to compute frozen contribution
            for (Size i=0; i<coll_idx; ++i) {
                Matrix3 Wmat =
                BpCollectionGeometry::frozen_jacobian_matrix_W(i, coll_idx,
                                                               bp_collection);
                Wmat.transpose();
                step_grads[i]._rotation += Wmat*frz_grad._translation;
            };

        };

    };


    // pencil function for electrostatics contributions to the gradients
    // those contributions are computed over all steps
    // the following code is from Juan Wei
    void add_electrostatics_contributions(StepGradientVec& step_grads,
                                          const BpCollection& bp_collection) {

        // base pairs
        const std::vector<BasePair>& bps = bp_collection.base_pairs();

        // sizing
        const Size n_steps = bp_collection.n_of_bp_steps();
        const Size n_bps = bp_collection.n_of_base_pairs();

        // loop over all steps
        for (Size k=0; k<n_steps; ++k) {

            // we go across the full set of base pairs
            for (Size i=0; i<k+1; ++i) {
                for (Size j=k+1; j<n_bps; ++j) {

                    if ((j-i)>(DNAElec::ExclusionRange-1)) {
                        const Vector3 Dij = bps[j].origin()-bps[i].origin();
                        const Real rij = Dij.norm();
                        step_grads[k]._translation +=
                        -DNAElec::dhConstant*std::exp(-DNAElec::kappaDebye*rij)*
                        (Real(1.0)/rij)*(DNAElec::kappaDebye+Real(1.0)/rij)*
                        Dij*(Real(1.0)/rij);
                    };

                };
            };

        };

    };


}


}

