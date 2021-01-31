// BpCollectionHessian functions
// Nicolas Clauvelin


#include "DNA/BpCollection.h"
#include "DNA/IndexManager.h"
#include "DNA/HessianFunctions.h"
#include "DNA/BpCollectionGeometry.h"
#include "DNA/BpCollectionHessian.h"


namespace BpCollectionHessian {


    // free collection Hessian matrix
    MatrixN free_collection_Hessian(const BpCollection& bp_collection,
                                    const IndexManager& idx_mgr) {

        // Hessian array
        const StepBlockArray2D h_mat =
        HessianFunctions::free_collection_Hessian_array(bp_collection, idx_mgr);

        // matrix
        return h_mat.regular_matrix();

    };


//// following code is in an unnamed namespace
//namespace {


    // EEDR boundary conditions B matrix as StepBlock
    StepBlock eedr_B_matrix(const Size& step_index,
                            const BpCollection& bp_collection,
                            const IndexManager& idx_mgr) {

        // bc reduced step index
        const Size& bc_idx = idx_mgr.bc_reduced_step().collection_index();

        // reduction matrix K
        const Matrix3 Kmat =
        Real(-1.0)*BpCollectionGeometry::jacobian_matrix_EER(step_index,
                                                             bc_idx,
                                                             bp_collection);

        // translation reduction matrix
        const Matrix3 mat = Real(-1.0)*Matrix3::identity_matrix();

        // step block assembly
        StepBlock block;
        block._rot_rot = Kmat;
        block._trans_trans = mat;

        return block;
        
    };


    // first order eedr contributions
    void order1_eerd_contributions(StepBlockArray2D& free_hessian,
                                   const BpCollection& bp_collection,
                                   const IndexManager& idx_mgr) {

        // sizing
        const Size n_steps = bp_collection.n_of_bp_steps();

        // bc reduced step collection index
        const Size& bc_idx = idx_mgr.bc_reduced_step().collection_index();

        // bc single step Hessian contribution
        const StepBlock& bc_h_block = free_hessian(bc_idx,bc_idx);

        // loop over all steps
        for (Size i=0; i<n_steps; ++i) {
            for (Size j=0; j<n_steps; ++j) {

                // exclude bc step
                if (i != bc_idx && j!= bc_idx) {

                    // eedr matrices
                    StepBlock BiT = eedr_B_matrix(i, bp_collection, idx_mgr);
                    BiT.transpose();
                    StepBlock Bj = eedr_B_matrix(j, bp_collection, idx_mgr);

                    // first correction
                    free_hessian(i,j) += BiT*bc_h_block*Bj;

                    // hessian blocks
                    const StepBlock& h_i_bc = free_hessian(i,bc_idx);
                    const StepBlock& h_j_bc = free_hessian(j,bc_idx);

                    // remaining corrections
                    free_hessian(i,j) += h_i_bc*Bj;
                    free_hessian(i,j) += BiT*h_j_bc;

                };

            };
        };

    };


//}

    
}
