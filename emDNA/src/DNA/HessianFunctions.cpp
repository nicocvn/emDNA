// HessianFunctions functions
// Nicolas Clauvelin


#include "DNA/IndexManager.h"
#include "DNA/BpCollection.h"
#include "DNA/BpCollectionGeometry.h"
#include "DNA/BpCollectionGradient.h"
#include "DNA/HessianFunctions.h"


namespace HessianFunctions {


    // function to symmetrized a force constant matrix
    StepBlock symmetrized_force_constants(const MatrixN& fmat) {

        MatrixN fmat_sym = Real(0.5)*fmat;
        fmat_sym.transpose();
        fmat_sym += Real(0.5)*fmat;

        // step block filling
        StepBlock step_fmat;
        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                step_fmat._rot_rot(i,j) = fmat_sym(i,j)*RAD_2_DEG*RAD_2_DEG;
        for (Size i=0; i<3; ++i)
            for (Size j=3; j<6; ++j)
                step_fmat._rot_trans(i,j-3) = fmat_sym(i,j)*RAD_2_DEG;
        for (Size i=3; i<6; ++i)
            for (Size j=0; j<3; ++j)
                step_fmat._trans_rot(i-3,j) = fmat_sym(i,j)*RAD_2_DEG;
        for (Size i=3; i<6; ++i)
            for (Size j=3; j<6; ++j)
                step_fmat._trans_trans(i-3,j-3) = fmat_sym(i,j);

        return step_fmat;

    };


    // free collection Hessian 2D map
    StepBlockArray2D
    free_collection_Hessian_array(const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr) {

        // data cache
        const BpDataCache dt_cache =
        bp_collection_data_cache(bp_collection, idx_mgr);

        // sizing
        const Size n_steps = bp_collection.n_of_bp_steps();

        // array
        StepBlockArray2D h_mat(n_steps);

        // diagonal terms
        for (Size i=0; i<n_steps; ++i) {

            // contrib
            const StepBlock diag_contrib =
            HessianFunctions::single_step_Hessian_contrib(i,
                                                          bp_collection,
                                                          idx_mgr,
                                                          dt_cache);

            // map insertion
            h_mat(i,i) = diag_contrib;

        };

        // upper off diagonal terms
        for (Size i=0; i<n_steps; ++i) {
            for (Size j=i+1; j<n_steps; ++j) {

                // contrib
                const StepBlock off_diag_contrib =
                HessianFunctions::coupling_Hessian_contrib(i, j,
                                                           bp_collection,
                                                           idx_mgr,
                                                           dt_cache);
                StepBlock off_diag_contrib_T = off_diag_contrib;
                off_diag_contrib_T.transpose();

                // map insertion
                h_mat(i,j) = off_diag_contrib;
                h_mat(j,i) = off_diag_contrib_T;
                
            };
        };
        
        return h_mat;

    };


    // Hessian diagonal contrib
    StepBlock single_step_Hessian_contrib(const Size& step_index,
                                          const BpCollection& bp_collection,
                                          const IndexManager& idx_mgr,
                                          const BpDataCache& dc) {

        // first order contrib
        StepBlock contrib =
        order1_diag_block_contrib(step_index, bp_collection, idx_mgr, dc);

        // PsiR contribution
        Matrix3 PsiR = order2_diag_block_contrib_PsiR(step_index,
                                                      bp_collection,
                                                      dc);
        contrib._rot_trans += PsiR;
        PsiR.transpose();
        contrib._trans_rot += PsiR;

        // PsiPsi contributions
        contrib._rot_rot +=
        order2_diag_block_contrib_PsiPsi(step_index, bp_collection, dc);
        contrib._rot_rot +=
        order2_diag_block_coupling_contrib_PsiPsi(step_index,
                                                  bp_collection, idx_mgr,
                                                  dc);

        return contrib;

    };


    // Hessian off diagonal coupling contrib
    StepBlock coupling_Hessian_contrib(const Size& step_i,
                                       const Size& step_j,
                                       const BpCollection& bp_collection,
                                       const IndexManager& idx_mgr,
                                       const BpDataCache& dc) {

        // ordering checking
        DS_ASSERT(step_i<step_j,
                  "computing Hessian first-order off diagonal term"
                  " with i>j");

        // first order contrib
        StepBlock contrib =
        order1_off_diag_block_contrib(step_i, step_j, bp_collection,
                                      dc);

        // PsiR contribution
        const Matrix3 PsiR =
        order2_off_diag_block_contrib_PsiR(step_i, step_j, bp_collection,
                                           dc);
        contrib._rot_trans += PsiR;

        // PsiPsi contribution
        const Matrix3 PsiPsi =
        order2_off_diag_block_contrib_PsiPsi(step_i, step_j,
                                             bp_collection,
                                             idx_mgr,
                                             dc);
        contrib._rot_rot += PsiPsi;

        return contrib;

    };


// following code is in an unnamed namespace
namespace {


    // symmetrized force constants matrix
    StepBlock symmetrized_force_constants(Size index,
                                          const BpCollection& bp_collection) {

        // force constants
        const MatrixN& fmat = bp_collection.bp_step_force_constants(index);
        MatrixN fmat_sym = Real(0.5)*(fmat+fmat.get_transpose());

        // step block filling
        StepBlock step_fmat;
        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                step_fmat._rot_rot(i,j) = fmat_sym(i,j)*RAD_2_DEG*RAD_2_DEG;
        for (Size i=0; i<3; ++i)
            for (Size j=3; j<6; ++j)
                step_fmat._rot_trans(i,j-3) = fmat_sym(i,j)*RAD_2_DEG;
        for (Size i=3; i<6; ++i)
            for (Size j=0; j<3; ++j)
                step_fmat._trans_rot(i-3,j) = fmat_sym(i,j)*RAD_2_DEG;
        for (Size i=3; i<6; ++i)
            for (Size j=3; j<6; ++j)
                step_fmat._trans_trans(i-3,j-3) = fmat_sym(i,j);

        return step_fmat;

    };


    // local jacobian matrix for single step
    StepBlock single_step_jacobian(Size step_index,
                                   const BpCollection& bp_collection) {

        const Matrix3& Rmat =
        BpCollectionGeometry::jacobian_matrix_R(step_index, bp_collection);
        const Matrix3& Tmat =
        BpCollectionGeometry::jacobian_matrix_T(step_index, bp_collection);

        return StepBlock(Matrix3::identity_matrix(), Rmat,
                         Matrix3(), Tmat);

    };


    // local jacobian matrix for coupled steps
    // index1 < index2
    StepBlock coupled_step_jacobian(Size index1, Size index2,
                                    const BpCollection& bp_collection) {

        const Matrix3& Umat =
        BpCollectionGeometry::jacobian_matrix_U(index1, index2, bp_collection);

        return StepBlock(Matrix3(), Umat,
                         Matrix3(), Matrix3());
        
    };











    // free collection first order contributions
    StepBlockArray2D
    free_collection_first_order(const BpCollection& bp_collection,
                                const IndexManager& idx_mgr) {

        // sizing
        const Size n_steps = bp_collection.n_of_bp_steps();

        // container
        // the container is initialized as zero
        StepBlockArray2D first_order_H(n_steps);

        // free step collection indices
        const std::vector<Size> free_idx =
        idx_mgr.free_step_collection_indices();
        const Size n_free_steps = free_idx.size();

        // loop over all free steps
        for (Size i=0; i<n_free_steps; ++i) {

            // collection index
            const Size& idx_i = free_idx[i];

            // step symmetrized force constants
            const StepBlock& fmat_i = symmetrized_force_constants(idx_i,
                                                                  bp_collection);

            // diagonal contribution - first part
            const StepBlock jmat_left =
            single_step_jacobian(idx_i, bp_collection);
            const StepBlock jmat_right = jmat_left.get_transpose();
            first_order_H(idx_i,idx_i) = jmat_left*fmat_i*jmat_right;

            // upper diagonal contributions
            for (Size j=i+1; j<n_free_steps; ++j) {

                // collection index
                const Size& idx_j = free_idx[j];

                // diagonal contribution - second part
                const StepBlock coupling_left =
                coupled_step_jacobian(idx_i, idx_j, bp_collection);
                const StepBlock coupling_right = coupling_left.get_transpose();
                first_order_H(idx_i,idx_i) += coupling_left*fmat_i*coupling_right;

//                // off-diagonal coupling contribution
//                const StepBlock
//                jjmat_right(Matrix3::identity_matrix(), Matrix3(),
//                            Rmat_jT, Tmat_jT);
//                hessian2D(idx_i,idx_j) = coupling_left*fsmat_i*jjmat_right;

            };

        };

        return first_order_H;

    };


    // bp collection data cache
    BpDataCache bp_collection_data_cache(const BpCollection& bp_collection,
                                         const IndexManager& idx_mgr) {

        // sizing
        const Size n_steps = bp_collection.n_of_bp_steps();

        // struct and memory allocation
        BpDataCache bp_coll_data;
        bp_coll_data._step_gradients.reserve(n_steps);
        bp_coll_data._R_matrices.reserve(n_steps);
        bp_coll_data._T_matrices.reserve(n_steps);
        bp_coll_data._S_vectors.reserve(n_steps);
        bp_coll_data._Lambda_vectors.reserve(n_steps);
        bp_coll_data._dLambda_matrices.reserve(n_steps);
        bp_coll_data._dXi_matrices.reserve(n_steps);

        // init
        for (Size i=0; i<n_steps; ++i) {

            // step gradient
            const BpStepParams& p =
            bp_collection.bp_step_params(i);
            const BpStepParams& p0 =
            bp_collection.bp_step_intrinsic_parameters(i);
            const MatrixN& fmat =
            bp_collection.bp_step_force_constants(i);
            bp_coll_data._step_gradients.
            push_back(BpCollectionGradient::
                      single_step_parameters_gradient(p, p0, fmat));

            // R and T jacobian matrices
            bp_coll_data._R_matrices.
            push_back(BpCollectionGeometry::
                      jacobian_matrix_R(i, bp_collection));
            bp_coll_data._T_matrices.
            push_back(BpCollectionGeometry::
                      jacobian_matrix_T(i, bp_collection));

            // S vectors
            const Vector3Triplet S_vectors =
            BpCollectionGeometry::rotation_vectors_S(i, bp_collection);
            bp_coll_data._S_vectors.
            push_back({S_vectors[0], S_vectors[1], S_vectors[2]});

            // euler angles
            const ZYZEulerAngles euler_angles =
            BpCollectionGeometry::
            ZYZEulerAngles_from_step_parameters(bp_collection.
                                                bp_step_params(i));

            // Lambda vectors
            const Vector3Triplet Lambda_vectors =
            BpGeometryFunctions::step_Lambda_vectors(euler_angles);
            bp_coll_data._Lambda_vectors.push_back(std::move(Lambda_vectors));

            // dLambda matrices
            const Matrix3Triplet dLambda_matrices =
            BpGeometryFunctions::step_dLambda_matrices(euler_angles);
            bp_coll_data._dLambda_matrices.
            push_back(std::move(dLambda_matrices));

            // dXi matrices
            const Matrix3Triplet dXi_matrices =
            BpGeometryFunctions::step_dXi_matrices(euler_angles);
            bp_coll_data._dXi_matrices.push_back(std::move(dXi_matrices));

        };
        
        return bp_coll_data;
        
    };


//    // free collection first order contributions function
//    // this function is used to create the step block 2D array
//    //
//    // this is the first order contribution diagonal and off diagonal
//    // corresponds to H_1^i and H_1^i,j in the notes
//    StepBlockArray2D
//    free_collection_first_order_contribs(const BpCollection& bp_collection,
//                                         const IndexManager& idx_mgr,
//                                         const BpDataCache& dc) {
//
//        // sizing
//        const Size n_steps = bp_collection.n_of_bp_steps();
//
//        // step block array
//        StepBlockArray2D hessian2D(n_steps);
//
//        // free steps collection indices
//        const std::vector<Size> free_idx =
//        idx_mgr.free_step_collection_indices();
//        const Size n_free_steps = free_idx.size();
//
//        // loop over all free steps
//        for (Size i=0; i<n_free_steps; ++i) {
//
//            // collection index
//            const Size& idx_i = free_idx[i];
//
//            // required data
//            const MatrixN& fmat_i =
//            bp_collection.bp_step_force_constants(idx_i);
//            const StepBlock& fsmat_i = symmetrized_force_constants(fmat_i);
//            const Matrix3& Rmat_i = dc._R_matrices[idx_i];
//            const Matrix3& Tmat_i = dc._T_matrices[idx_i];
//
//            // diagonal contribution - first part
//            const StepBlock jmat_left(Matrix3::identity_matrix(), Rmat_i,
//                                      Matrix3(), Tmat_i);
//            const StepBlock jmat_right = jmat_left.get_transpose();
//            hessian2D(idx_i,idx_i) = jmat_left*fsmat_i*jmat_right;
//
//            // upper diagonal contributions
//            for (Size j=i+1; j<n_free_steps; ++j) {
//
//                // collection index
//                const Size& idx_j = free_idx[j];
//
//                // required data
//                const Matrix3 Rmat_jT = dc._R_matrices[idx_j].get_transpose();
//                const Matrix3 Tmat_jT = dc._T_matrices[idx_j].get_transpose();
//                const Matrix3 Umat_ij =
//                BpCollectionGeometry::jacobian_matrix_U(idx_i, idx_j,
//                                                        bp_collection);
//
//                // diagonal contribution - second part
//                const StepBlock coupling_left(Matrix3(), Umat_ij,
//                                              Matrix3(), Matrix3());
//                const StepBlock coupling_right = coupling_left.get_transpose();
//                hessian2D(idx_i,idx_i) += coupling_left*fsmat_i*coupling_right;
//
//                // coupling contribution
//                const StepBlock
//                jjmat_right(Matrix3::identity_matrix(), Matrix3(),
//                            Rmat_jT, Tmat_jT);
//                hessian2D(idx_i,idx_j) = coupling_left*fsmat_i*jjmat_right;
//
//            };
//
//        };
//
//        return hessian2D;
//
//    };
//
//
//    // free collection second order contributions function
//    void
//    free_collection_second_order_contribs(StepBlockArray2D& hessian2D,
//                                          const BpCollection& bp_collection,
//                                          const IndexManager& idx_mgr,
//                                          const BpDataCache& dc) {
//
//        // sizing
//        const Size n_steps = bp_collection.n_of_bp_steps();
//
//        // free steps collection indices
//        const std::vector<Size> free_idx =
//        idx_mgr.free_step_collection_indices();
//        const Size n_free_steps = free_idx.size();
//
//        // loop over all free steps
//        for (Size i=0; i<n_free_steps; ++i) {
//
//            // collection index
//            const Size& idx_i = free_idx[i];
//
//            // required data
//            const BpStepParams& p_i = bp_collection.bp_step_params(idx_i);
//            const Vector3 rho_vec_i = {
//                p_i.value(SHIFT), p_i.value(SLIDE), p_i.value(RISE)
//            };
//            const StepGradient& grad_i = dc._step_gradients[idx_i];
//            const Matrix3& Rmat_i = dc._R_matrices[idx_i];
//            const Matrix3& Tmat_i = dc._T_matrices[idx_i];
//            const Vector3Triplet& Lambda_i = dc._Lambda_vectors[idx_i];
//            const Matrix3Triplet dLambda_i = {
//                dc._dLambda_matrices[idx_i][0].get_transpose(),
//                dc._dLambda_matrices[idx_i][1].get_transpose(),
//                dc._dLambda_matrices[idx_i][2].get_transpose()
//            };
//
//            // diagonal contribution - psi/r
//            hessian2D(idx_i,idx_i)._rot_trans +=
//            order2_diag_block_PsiR_contrib(grad_i, Tmat_i, Lambda_i);
//
//            // diagonal contribution - psi/psi
//            hessian2D(idx_i,idx_i)._rot_rot +=
//            order2_diag_block_PsiPsi_contrib1(grad_i, Rmat_i, Lambda_i);
//            hessian2D(idx_i,idx_i)._rot_rot +=
//            order2_diag_block_PsiPsi_contrib2(grad_i, rho_vec_i, dLambda_i);
//
//            // upper diagonal contributions
//            for (Size j=i+1; j<n_free_steps; ++j) {
//
//                // collection index
//                const Size& idx_j = free_idx[j];
//
//                // required data
//                const StepGradient& grad_j = dc._step_gradients[idx_j];
//                const Matrix3 Rmat_jT = dc._R_matrices[idx_j].get_transpose();
//                const Matrix3 Tmat_jT = dc._T_matrices[idx_j].get_transpose();
//                const Matrix3 Umat_ij =
//                BpCollectionGeometry::jacobian_matrix_U(idx_i, idx_j,
//                                                        bp_collection);
//                const Vector3Triplet& Lambda_j = dc._Lambda_vectors[idx_j];
//
//                
//                
//            };
//            
//        };
//
//    };


    // first-order Hessian diagonal contribution
    StepBlock order1_diag_block_contrib(const Size& step_index,
                                        const BpCollection& bp_collection,
                                        const IndexManager& idx_mgr,
                                        const BpDataCache& dc) {

        // force constant matrix
        const StepBlock fmat =
        symmetrized_force_constants(step_index, bp_collection);

        // local jacobian matrices
        const Matrix3& Rmat = dc._R_matrices[step_index];
        const Matrix3& Tmat = dc._T_matrices[step_index];

        // local transformation matrices
        StepBlock trans_mat_left(Matrix3::identity_matrix(), Rmat,
                                 Matrix3(), Tmat);
        StepBlock trans_mat_right = trans_mat_left;
        trans_mat_right.transpose();

        // local contribution
        const StepBlock local_contrib = trans_mat_left*fmat*trans_mat_right;

        // free steps iterator bounds
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // loop over downward free steps for coupling contributions
        // frozen steps do not contribute
        StepBlock coupling_contrib;
        for (auto it = free_begin; it != free_end; ++it) {

            // collection index
            const Size& coll_idx = it->collection_index();

            // check if it is a downward step
            if (coll_idx > step_index) {

                // jacobian matrix
                const Matrix3 Umat =
                BpCollectionGeometry::jacobian_matrix_U(step_index,
                                                        coll_idx,
                                                        bp_collection);

                // transformation matrix
                StepBlock coupling_left(Matrix3(), Umat,
                                        Matrix3(), Matrix3());
                StepBlock coupling_right = coupling_left;
                coupling_right.transpose();
                
                // contribution
                coupling_contrib += coupling_left*fmat*coupling_right;
                
            };
            
        };
        
        return local_contrib+coupling_contrib;

    };


    // second-order Hessian diagonal contribution - rot/trans
    Matrix3 order2_diag_block_PsiR_contrib(const StepGradient& grad,
                                           const Matrix3& Tmat,
                                           const Vector3Triplet& Lambda) {

        // matrix assembly
        Matrix3 A_PsiR;
        A_PsiR.set_row(0, Tmat*(grad.
                                _translation.cross(Lambda[0])));
        A_PsiR.set_row(1, Tmat*(grad.
                                _translation.cross(Lambda[1])));
        A_PsiR.set_row(2, Tmat*(grad.
                                _translation.cross(Lambda[2])));

        return A_PsiR;

    };
    Matrix3 order2_diag_block_contrib_PsiR(const Size& step_index,
                                           const BpCollection& bp_collection,
                                           const BpDataCache& dc) {

        // jacobian matrix T
        const Matrix3& Tmat = dc._T_matrices[step_index];

        // Lambda vectors
        const Vector3Triplet& Lambda_vectors = dc._Lambda_vectors[step_index];

        // step parameters gradient
        const StepGradient& step_gradient = dc._step_gradients[step_index];

        // matrix assembly
        Matrix3 A_PsiR;
        A_PsiR.set_row(0, Tmat*(step_gradient.
                                _translation.cross(Lambda_vectors[0])));
        A_PsiR.set_row(1, Tmat*(step_gradient.
                                _translation.cross(Lambda_vectors[1])));
        A_PsiR.set_row(2, Tmat*(step_gradient.
                                _translation.cross(Lambda_vectors[2])));

        return A_PsiR;

    };


    // second-order Hessian diagonal contribution - rot/rot
    Matrix3 order2_diag_block_PsiPsi_contrib1(const StepGradient& grad,
                                              const Matrix3& Rmat,
                                              const Vector3Triplet& Lambda) {

        // matrix assembly
        Matrix3 A_PsiPsi_1;
        A_PsiPsi_1.set_row(0, Rmat*(grad.
                                    _translation.cross(Lambda[0])));
        A_PsiPsi_1.set_row(1, Rmat*(grad.
                                    _translation.cross(Lambda[1])));
        A_PsiPsi_1.set_row(2, Rmat*(grad.
                                    _translation.cross(Lambda[2])));

        return A_PsiPsi_1;

    };
    Matrix3 order2_diag_block_PsiPsi_contrib2(const StepGradient& grad,
                                              const Vector3& rho_vec,
                                              const Matrix3Triplet& dLambda) {

        // matrix assembly
        Matrix3 A_PsiPsi_2;
        A_PsiPsi_2.set_row(0,
                           dLambda[0]*(rho_vec.
                                       cross(grad._translation)));
        A_PsiPsi_2.set_row(1,
                           dLambda[1]*(rho_vec.
                                       cross(grad._translation)));
        A_PsiPsi_2.set_row(2,
                           dLambda[2]*(rho_vec.
                                       cross(grad._translation)));

        return A_PsiPsi_2;

    };
    Matrix3
    order2_diag_block_contrib_PsiPsi(const Size& step_index,
                                     const BpCollection& bp_collection,
                                     const BpDataCache& dc) {

        // jacobian matrix R
        const Matrix3& Rmat = dc._R_matrices[step_index];

        // Lambda vectors
        const Vector3Triplet& Lambda_vectors = dc._Lambda_vectors[step_index];

        // Lambda vectors derivatives
        Matrix3Triplet dLambda_matsT = dc._dLambda_matrices[step_index];
        for (Size i=0; i<3; ++i)
            dLambda_matsT[i].transpose();

        // step parameters gradient
        const BpStepParams p = bp_collection.bp_step_params(step_index);
        const StepGradient& step_gradient = dc._step_gradients[step_index];

        // rho parameters vector
        const Vector3 rho_vec(p.value(SHIFT), p.value(SLIDE), p.value(RISE));

        // first part
        Matrix3 A_PsiPsi_1;
        A_PsiPsi_1.set_row(0, Rmat*(step_gradient.
                                    _translation.cross(Lambda_vectors[0])));
        A_PsiPsi_1.set_row(1, Rmat*(step_gradient.
                                    _translation.cross(Lambda_vectors[1])));
        A_PsiPsi_1.set_row(2, Rmat*(step_gradient.
                                    _translation.cross(Lambda_vectors[2])));

        // second part
        Matrix3 A_PsiPsi_2;
        A_PsiPsi_2.set_row(0,
                           dLambda_matsT[0]*(rho_vec.
                                            cross(step_gradient._translation)));
        A_PsiPsi_2.set_row(1,
                           dLambda_matsT[1]*(rho_vec.
                                            cross(step_gradient._translation)));
        A_PsiPsi_2.set_row(2,
                           dLambda_matsT[2]*(rho_vec.
                                            cross(step_gradient._translation)));

        return A_PsiPsi_1+A_PsiPsi_2;

    };


    // second-order Hessian coupling diagonal contribution
    Matrix3
    order2_diag_block_coupling_contrib_PsiPsi(const Size& step_index,
                                              const BpCollection& bp_collection,
                                              const IndexManager& idx_mgr,
                                              const BpDataCache& dc) {

        // step frame matrix
        const Matrix3 frame_mat =
        bp_collection.base_pair(step_index).column_axes();

        // Xi matrix derivatives
        const Matrix3Triplet& dXi_mats = dc._dXi_matrices[step_index];

        // local S vectors
        const Vector3Triplet& S_vectors_i = dc._S_vectors[step_index];

        // free steps iterator bounds
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // loop over free steps
        Matrix3 Ap_PsiPsi_1, Ap_PsiPsi_2;
        for (auto it = free_begin; it != free_end; ++it) {

            // collection index
            const Size& coll_idx = it->collection_index();

            // downward step
            if (coll_idx > step_index) {

                // jacobian matrix T
                Matrix3 Tmat_j = dc._T_matrices[coll_idx];
                Tmat_j.transpose();

                // transformation matrices
                std::vector<Matrix3> tr_mats;
                for (Size i=0; i<3; ++i) {
                    tr_mats.push_back(Tmat_j*frame_mat*dXi_mats[i]);
                    tr_mats.back().transpose();
                };

                // step parameters gradient
                const BpStepParams& p = bp_collection.bp_step_params(coll_idx);
                const StepGradient& step_gradient =
                dc._step_gradients[coll_idx];

                // rho parameters vector
                const Vector3
                rho_vec(p.value(SHIFT), p.value(SLIDE), p.value(RISE));

                // first part
                Matrix3 Ap_PsiPsi_1_add;
                Ap_PsiPsi_1_add.set_col(0,
                                    tr_mats[0]*
                                    (step_gradient.
                                     _translation.cross(rho_vec)));
                Ap_PsiPsi_1_add.set_col(1,
                                    tr_mats[1]*
                                    (step_gradient.
                                     _translation.cross(rho_vec)));
                Ap_PsiPsi_1_add.set_col(2,
                                    tr_mats[2]*
                                    (step_gradient.
                                     _translation.cross(rho_vec)));

                // transformed S vectors
                std::vector<Vector3> trans_S_vectors;
                for (Size i=0; i<3; ++i) {
                    trans_S_vectors.push_back(Tmat_j*S_vectors_i[i]);
                };

                // second part
                Matrix3 Ap_PsiPsi_2_add;
                for (Size i=0; i<3; ++i) {
                    for (Size j=0; j<3; ++j) {
                        const Vector3 v =
                        (rho_vec.
                         cross(trans_S_vectors[i])).cross(trans_S_vectors[j]);
                        Ap_PsiPsi_2_add(i,j) = v.dot(step_gradient.
                                                     _translation);
                    };
                };

                // accumulation
                Ap_PsiPsi_1 += Ap_PsiPsi_1_add;
                Ap_PsiPsi_2 += Ap_PsiPsi_2_add;

            };

        };

        return Ap_PsiPsi_1 + Ap_PsiPsi_2;

    };


    // first-order Hessian off-diagonal contribution
    StepBlock order1_off_diag_block_contrib(const Size& step_i,
                                            const Size& step_j,
                                            const BpCollection& bp_collection,
                                            const BpDataCache& dc) {

        // force constant matrix for step i
        const StepBlock fmat_i =
        symmetrized_force_constants(step_i, bp_collection);

        // local matrices
        Matrix3 Rmat_j = dc._R_matrices[step_j];
        Matrix3 Tmat_j = dc._T_matrices[step_j];
        Rmat_j.transpose();
        Tmat_j.transpose();
        const Matrix3 Umat_ij =
        BpCollectionGeometry::jacobian_matrix_U(step_i, step_j, bp_collection);

        // transformation matrices
        const StepBlock trans_mat_left(Matrix3(), Umat_ij,
                                       Matrix3(), Matrix3());

        const StepBlock trans_mat_right(Matrix3::identity_matrix(), Matrix3(),
                                        Rmat_j, Tmat_j);

        // Lambda vectors - j
        const Vector3Triplet& Lambda_vecs = dc._Lambda_vectors[step_j];

        // step parameters gradient
        const StepGradient& step_gradient = dc._step_gradients[step_j];

        // matrix assembly
        Matrix3 PsiPsi;
        PsiPsi.set_col(0,
                       Umat_ij*(step_gradient.
                                _translation.cross(Lambda_vecs[0])));
        PsiPsi.set_col(1,
                       Umat_ij*(step_gradient.
                                _translation.cross(Lambda_vecs[1])));
        PsiPsi.set_col(2,
                       Umat_ij*(step_gradient.
                                _translation.cross(Lambda_vecs[2])));

        // first order contrib
        StepBlock contrib = trans_mat_left*fmat_i*trans_mat_right;
        contrib._rot_rot += PsiPsi;

        return contrib;

    };


    // second-order Hessian off-diagonal contribution - rot/trans
    Matrix3 order2_off_diag_block_PsiR_contrib(const StepGradient& grad_j,
                                               const Matrix3& Tmat_j,
                                               const Vector3Triplet& Svecs_i) {

        // matrix assembly
        Matrix3 A_ij_PsiR;
        A_ij_PsiR.set_row(0,
                          Svecs_i[0].
                          cross(Tmat_j*grad_j._translation));
        A_ij_PsiR.set_row(1,
                          Svecs_i[1].
                          cross(Tmat_j*grad_j._translation));
        A_ij_PsiR.set_row(2,
                          Svecs_i[2].
                          cross(Tmat_j*grad_j._translation));

        return A_ij_PsiR;

    };
    Matrix3
    order2_off_diag_block_contrib_PsiR(const Size& step_i, const Size& step_j,
                                       const BpCollection& bp_collection,
                                       const BpDataCache& dc) {

        // S vectors - i
        const Vector3Triplet& Svectors = dc._S_vectors[step_i];

        // T matrix - j
        const Matrix3& Tmat = dc._T_matrices[step_j];

        // step parameters gradient - j
        const StepGradient step_gradient = dc._step_gradients[step_j];

        // matrix assembly
        Matrix3 A_ij_PsiR;
        A_ij_PsiR.set_row(0,
                          Svectors[0].
                          cross(Tmat*step_gradient._translation));
        A_ij_PsiR.set_row(1,
                          Svectors[1].
                          cross(Tmat*step_gradient._translation));
        A_ij_PsiR.set_row(2,
                          Svectors[2].
                          cross(Tmat*step_gradient._translation));

        return A_ij_PsiR;

    };


    // second-order-order Hessian off-diagonal contribution - rot/rot
    Matrix3
    order2_off_diag_block_contrib_PsiPsi(const Size& step_i, const Size& step_j,
                                         const BpCollection& bp_collection,
                                         const IndexManager& idx_mgr,
                                         const BpDataCache& dc) {

        // force constant matrix for step i
        const StepBlock fmat_i =
        symmetrized_force_constants(step_i, bp_collection);

        // S vectors - i
        const Vector3Triplet& Svectors_i = dc._S_vectors[step_i];

        // S vectors - j
        const Vector3Triplet& Svectors_j = dc._S_vectors[step_j];

        // free steps iterator bounds
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // loop on downward free steps
        Matrix3 A_ij_PsiPsi1_acc, A_ij_PsiPsi2_acc;
        for (auto it = free_begin; it != free_end; ++ it) {

            // collection index
            const Size& coll_idx = it->collection_index();

            // downward step
            if (coll_idx > step_j) {

                // step parameters
                const BpStepParams& p_h =
                bp_collection.bp_step_params(coll_idx);
                const Vector3 rho_vec_h = {
                    p_h.value(SHIFT), p_h.value(SLIDE), p_h.value(RISE)
                };

                // step parameters gradient
                const StepGradient& grad_h = dc._step_gradients[coll_idx];

                // jacobian matrix T
                Matrix3 Tmat_h = dc._T_matrices[coll_idx];
                Tmat_h.transpose();

                // jacobian matrices U
                const Matrix3 U_ih =
                BpCollectionGeometry::jacobian_matrix_U(step_i, coll_idx,
                                                        bp_collection);
                Matrix3 U_jh =
                BpCollectionGeometry::jacobian_matrix_U(step_j, coll_idx,
                                                        bp_collection);
                U_jh.transpose();

                // matrix assembly - 1
                Matrix3 m1 = U_ih*(fmat_i._trans_trans)*U_jh;

                // matrix assembly - 2
                Matrix3 m2;
                for (Size k=0; k<3; ++k) {
                    for (Size l=0; l<3; ++l) {
                        const Vector3 ts1 = Tmat_h*Svectors_j[l];
                        const Vector3 ts2 = Tmat_h*Svectors_i[k];
                        m2(k,l) =
                        (ts1.cross(ts2.cross(rho_vec_h))).
                        dot(grad_h._translation);
                    };
                };

                // accumulation
                A_ij_PsiPsi1_acc += m1;
                A_ij_PsiPsi2_acc += m2;

            };

        };

        return A_ij_PsiPsi1_acc+A_ij_PsiPsi2_acc;

    };


}


}
