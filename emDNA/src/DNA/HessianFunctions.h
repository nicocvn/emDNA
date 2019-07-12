// HessianFunctions functions
// Nicolas Clauvelin


// functions for the computation of the Hessian
//
// there is a lot of room for optimization here ...


#ifndef emDNA_HessianFunctions_h
#define emDNA_HessianFunctions_h


#include <StepBlock.h>
class StepGradient;
class BpCollection;
class IndexManager;
class BpCacheData;


namespace HessianFunctions {


    namespace {
        // data structure for caching
        using Vector3Triplet = std::vector<Vector3>;
        using Matrix3Triplet = std::vector<Matrix3>;
        struct BpDataCache {
            std::vector<StepGradient> _step_gradients;
            std::vector<Matrix3> _R_matrices;
            std::vector<Matrix3> _T_matrices;
            std::vector<Vector3Triplet> _S_vectors;
            std::vector<Vector3Triplet> _Lambda_vectors;
            std::vector<Matrix3Triplet> _dLambda_matrices;
            std::vector<Matrix3Triplet> _dXi_matrices;
        };
    }

    // function to symmetrized a force constant matrix
    StepBlock symmetrized_force_constants(const MatrixN& fmat);

    // free collection Hessian 2D array
    StepBlockArray2D
    free_collection_Hessian_array(const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr);

    // Hessian diagonal contrib
    StepBlock single_step_Hessian_contrib(const Size& step_index,
                                          const BpCollection& bp_collection,
                                          const IndexManager& idx_mgr,
                                          const BpDataCache& dc);

    // Hessian off diagonal coupling contrib
    StepBlock coupling_Hessian_contrib(const Size& step_i,
                                       const Size& step_j,
                                       const BpCollection& bp_collection,
                                       const IndexManager& idx_mgr,
                                       const BpDataCache& dc);
    

    // following code is in an unnamed namespace
    namespace {


    // symmetrized force constants matrix
    StepBlock symmetrized_force_constants(Size index,
                                          const BpCollection& bp_collection);

    // local jacobian matrix for single step
    StepBlock single_step_jacobian(Size step_index,
                                   const BpCollection& bp_collection);

    // local jacobian matrix for coupled steps
    // index1 < index2
    StepBlock coupled_step_jacobian(Size index1, Size index2,
                                    const BpCollection& bp_collection);






        


    // free collection first order contributions
    StepBlockArray2D
    free_collection_first_order(const BpCollection& bp_collection,
                                const IndexManager& idx_mgr);


    // bp collection data cache
    BpDataCache bp_collection_data_cache(const BpCollection& bp_collection,
                                         const IndexManager& idx_mgr);

    // free collection first order contributions function
    // this function is used to create the step block 2D array
    StepBlockArray2D
    free_collection_first_order_contribs(const BpCollection& bp_collection,
                                         const IndexManager& idx_mgr,
                                         const BpDataCache& dc);

    // free collection second order contributions function
    void
    free_collection_second_order_contribs(StepBlockArray2D& hessian2D,
                                          const BpCollection& bp_collection,
                                          const IndexManager& idx_mgr,
                                          const BpDataCache& dc);


    StepBlock order1_diag_block_contrib(const Size& step_index,
                                        const BpCollection& bp_collection,
                                        const IndexManager& idx_mgr,
                                        const BpDataCache& dc);

    // second-order Hessian diagonal contribution - rot/trans
    Matrix3 order2_diag_block_PsiR_contrib(const StepGradient& grad,
                                           const Matrix3& Tmat,
                                           const Vector3Triplet& Lambda);
    Matrix3 order2_diag_block_contrib_PsiR(const Size& step_index,
                                           const BpCollection& bp_collection,
                                           const BpDataCache& dc);

    // second-order Hessian diagonal contribution - rot/rot
    Matrix3 order2_diag_block_PsiPsi_contrib1(const StepGradient& grad,
                                              const Matrix3& Rmat,
                                              const Vector3Triplet& Lambda);
    Matrix3 order2_diag_block_PsiPsi_contrib2(const StepGradient& grad,
                                              const Vector3& rho_vec,
                                              const Matrix3Triplet& dLambda);
    Matrix3
    order2_diag_block_contrib_PsiPsi(const Size& step_index,
                                     const BpCollection& bp_collection,
                                     const BpDataCache& dc);

    // second-order Hessian coupling diagonal contribution
    Matrix3
    order2_diag_block_coupling_contrib_PsiPsi(const Size& step_index,
                                              const BpCollection& bp_collection,
                                              const IndexManager& idx_mgr,
                                              const BpDataCache& dc);

    // first-order Hessian off-diagonal contribution
    StepBlock order1_off_diag_block_contrib(const Size& step_i,
                                            const Size& step_j,
                                            const BpCollection& bp_collection,
                                            const BpDataCache& dc);

    // second-order Hessian off-diagonal contribution - rot/trans
    Matrix3 order2_off_diag_block_PsiR_contrib(const StepGradient& grad_j,
                                               const Matrix3& Tmat_j,
                                               const Vector3Triplet& Svecs_i);
    Matrix3
    order2_off_diag_block_contrib_PsiR(const Size& step_i, const Size& step_j,
                                       const BpCollection& bp_collection,
                                       const BpDataCache& dc);

    // second-order-order Hessian off-diagonal contribution - rot/rot
    Matrix3
    order2_off_diag_block_contrib_PsiPsi(const Size& step_i, const Size& step_j,
                                         const BpCollection& bp_collection,
                                         const IndexManager& idx_mgr,
                                         const BpDataCache& dc);


    }


}


#endif  // emDNA_HessianFunctions_h
