// BpCollectionGeometry functions
// Nicolas Clauvelin


#include <StepBlock.h>
#include <BpCollectionGeometry.h>


namespace BpCollectionGeometry {


    // step parameters to ZYZ Euler angles
    ZYZEulerAngles ZYZEulerAngles_from_step_parameters(const BpStepParams& p) {
        return Theta_2_ZYZEuler(ThetaAngles(p.value(TILT),
                                            p.value(ROLL),
                                            p.value(TWIST)));
    };


    // ZYZ Euler angles from step dofs
    ZYZEulerAngles ZYZEulerAngles_from_step_dofs(const BpStepDofs& dofs) {
        return Theta_2_ZYZEuler(ThetaAngles(dofs.value(TILTrad)*RAD_2_DEG,
                                            dofs.value(ROLLrad)*RAD_2_DEG,
                                            dofs.value(TWISTrad)*RAD_2_DEG));
    };


    // step frame
    Matrix3 step_frame_column_axes(Size step_index,
                                   const BpCollection& bp_collection) {

        // step ZYZ Euler angles
        const ZYZEulerAngles euler_angles =
        ZYZEulerAngles_from_step_parameters(bp_collection.
                                            bp_step_params(step_index));

        // step frame column axes
        return Matrix3(bp_collection.base_pair(step_index).column_axes()*
                       step_frame_rotation_matrix(euler_angles));

    };


    // jacobian matrix T
    // this corresponds to the transpose of the matrix R defined in the notes
    Matrix3 jacobian_matrix_T(Size step_index,
                              const BpCollection& bp_collection) {
        return step_frame_column_axes(step_index, bp_collection);
    };


    // jacobian matrix R
    // this corresponds to the transpose of the matrix R defined in the notes
    Matrix3 jacobian_matrix_R(Size step_index,
                              const BpCollection& bp_collection) {

        // step ZYZ Euler angles
        const ZYZEulerAngles euler_angles =
        ZYZEulerAngles_from_step_parameters(bp_collection.
                                            bp_step_params(step_index));

        // Lambda vectors
        const std::vector<Vector3> Lambda_vectors =
        step_Lambda_vectors(euler_angles);

        // rho vector
        const Vector3 rhovec =
        Vector3(bp_collection.bp_step_params(step_index).value(SHIFT),
                bp_collection.bp_step_params(step_index).value(SLIDE),
                bp_collection.bp_step_params(step_index).value(RISE));

        // matrix components
        const Vector3 v1(Lambda_vectors[0].cross(rhovec));
        const Vector3 v2(Lambda_vectors[1].cross(rhovec));
        const Vector3 v3(Lambda_vectors[2].cross(rhovec));

        // matrix R
        Matrix3 Rmat;
        Rmat.set_row(0, v1);
        Rmat.set_row(1, v2);
        Rmat.set_row(2, v3);

        return Rmat;
        
    };


    // infinitesimal rotation vectors S
    std::vector<Vector3> rotation_vectors_S(Size step_index,
                                            const BpCollection& bp_collection) {

        // step ZYZ Euler angles
        const ZYZEulerAngles euler_angles =
        ZYZEulerAngles_from_step_parameters(bp_collection.
                                            bp_step_params(step_index));

        // Xi matrix
        const Matrix3 Ximat(BpGeometryFunctions::step_Xi_matrix(euler_angles));

        // bp frame
        const Matrix3 dmat(bp_collection.base_pair(step_index).column_axes());

        // packing
        std::vector<Vector3> S_vectors;
        S_vectors.reserve(3);
        S_vectors.push_back(dmat*Ximat.col(0));
        S_vectors.push_back(dmat*Ximat.col(1));
        S_vectors.push_back(dmat*Ximat.col(2));

        return S_vectors;

    };


    // jacobian matrix U
    // index1 < index2
    // this corresponds to the transpose of the matrix U defined in the notes
    Matrix3 jacobian_matrix_U(Size index1, Size index2,
                              const BpCollection& bp_collection) {

        DS_ASSERT(index1<index2,
                  "computing U(i,j) with i>j");

        // rho vector
        const Vector3 rhovec(bp_collection.bp_step_params(index2).value(SHIFT),
                             bp_collection.bp_step_params(index2).value(SLIDE),
                             bp_collection.bp_step_params(index2).value(RISE));

        // S vectors
        const std::vector<Vector3> S_vectors(rotation_vectors_S(index1,
                                                                bp_collection));

        // step frame column axes
        Matrix3 step_frame_mat_T = step_frame_column_axes(index2,
                                                          bp_collection);
        step_frame_mat_T.transpose();

        // matrix components
        const Vector3 v1(rhovec.cross(step_frame_mat_T*S_vectors[0]));
        const Vector3 v2(rhovec.cross(step_frame_mat_T*S_vectors[1]));
        const Vector3 v3(rhovec.cross(step_frame_mat_T*S_vectors[2]));

        // matrix U
        Matrix3 Umat;
        Umat.set_row(0, v1);
        Umat.set_row(1, v2);
        Umat.set_row(2, v3);

        return Umat;
        
    };


    // jacobian matrix J
    // J^index1,index2
    // index1 < index2 (index2 has to be related to the bc reduced step)
    Matrix3 jacobian_matrix_rotS(Size index1, Size index2,
                                 const BpCollection& bp_collection) {

        DS_ASSERT(index1<index2,
                  "computing rotS(i,j) with i>j");

        // S vectors
        const std::vector<Vector3> S_vectors =
        rotation_vectors_S(index1, bp_collection);
        Matrix3 Smat;
        Smat.set_col(0, S_vectors[0]);
        Smat.set_col(1, S_vectors[1]);
        Smat.set_col(2, S_vectors[2]);

        // bp frame column axes matrices
        Matrix3 dmatT(bp_collection.base_pair(index2).column_axes());
        dmatT.transpose();

        return dmatT*Smat;

    };


    // jacobian matrix SigmaJ
    // (Omega^index1).J^(index1)_(index2)
    // index1 < index2 (index2 has to be related to the bc reduced step)
    Matrix3 jacobian_matrix_EER(Size index1, Size index2,
                                   const BpCollection& bp_collection) {

        DS_ASSERT(index1<index2,
                  "computing EER(i,j)=K with i>j");

        // step ZYZ Euler angles
        const ZYZEulerAngles euler_angles =
        ZYZEulerAngles_from_step_parameters(bp_collection.
                                            bp_step_params(index2));

        // matrix Sigma
        Matrix3 Omega(step_Omega_matrix(euler_angles));

        return Omega*jacobian_matrix_rotS(index1, index2, bp_collection);

    };


    // frozen step jacobian matrices W and C
    // index1 < index2
    Matrix3 frozen_jacobian_matrix_W(Size index1, Size index2,
                                     const BpCollection& bp_collection) {

        DS_ASSERT(index1<index2,
                  "computing W(i,j) with i>j");


        // r dofs for index2
        const Vector3 r(bp_collection.bp_step_dofs(index2).value(R1),
                        bp_collection.bp_step_dofs(index2).value(R2),
                        bp_collection.bp_step_dofs(index2).value(R3));

        // S vectors
        const std::vector<Vector3> S_vectors(rotation_vectors_S(index1,
                                                                bp_collection));

        // W matrix components
        const Vector3 v1(S_vectors[0].cross(r));
        const Vector3 v2(S_vectors[1].cross(r));
        const Vector3 v3(S_vectors[2].cross(r));

        // W matrix
        Matrix3 Wmat;
        Wmat.set_col(0, v1);
        Wmat.set_col(1, v2);
        Wmat.set_col(2, v3);

        return Wmat;

    };


}
