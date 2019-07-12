// BpCollectionGeometry functions
// Nicolas Clauvelin


// functions for the geometry of bp collections
// see implementation notes for details


#ifndef emDNA_BpCollectionGeometry_h
#define emDNA_BpCollectionGeometry_h


#include <BpGeometryFunctions.h>
#include <BpCollection.h>
using namespace BpGeometryFunctions;
class StepBlock;


namespace BpCollectionGeometry {
    

    // ZYZ Euler angles from step parameters
    ZYZEulerAngles ZYZEulerAngles_from_step_parameters(const BpStepParams& p);

    // ZYZ Euler angles from step dofs
    ZYZEulerAngles ZYZEulerAngles_from_step_dofs(const BpStepDofs& dofs);

    // step frame
    Matrix3 step_frame_column_axes(Size step_index,
                                   const BpCollection& bp_collection);

    // jacobian matrix T
    // this corresponds to the transpose of the matrix T defined in the notes
    Matrix3 jacobian_matrix_T(Size step_index,
                              const BpCollection& bp_collection);

    // jacobian matrix R
    // this corresponds to the transpose of the matrix R defined in the notes
    Matrix3 jacobian_matrix_R(Size step_index,
                              const BpCollection& bp_collection);

    // infinitesimal rotation vectors S
    std::vector<Vector3> rotation_vectors_S(Size step_index,
                                            const BpCollection& bp_collection);

    // jacobian matrix U
    // index1 < index2
    // this corresponds to the transpose of the matrix U defined in the notes
    Matrix3 jacobian_matrix_U(Size index1, Size index2,
                              const BpCollection& bp_collection);

    // jacobian matrix rotS
    // J^index1_index2
    // index1 < index2 (index2 has to be related to the bc reduced step)
    // this corresponds to a K without Omega
    Matrix3 jacobian_matrix_rotS(Size index1, Size index2,
                                 const BpCollection& bp_collection);

    // jacobian matrix EER
    // (Omega^index1).J^(index1)_(index2)
    // index1 < index2 (index2 has to be related to the bc reduced step)
    // this corresponds to the matrix K in the notes
    Matrix3 jacobian_matrix_EER(Size index1, Size index2,
                                   const BpCollection& bp_collection);

    // frozen step jacobian matrix W
    // index1 < index2
    Matrix3 frozen_jacobian_matrix_W(Size index1, Size index2,
                                     const BpCollection& bp_collection);


}


#endif  // emDNA_BpCollectionGeometry_h
