// BpGeometryFunctions functions
// Nicolas Clauvelin


// functions for the geometry of bp steps


#ifndef emDNA_BpGeometryFunctions_h
#define emDNA_BpGeometryFunctions_h


#include <emDNA_Includes.h>


namespace BpGeometryFunctions {


    // ZYZ Euler angles structure
    struct ZYZEulerAngles {
        Real _ZetaRadian;
        Real _KappaRadian;
        Real _EtaRadian;
        ZYZEulerAngles() :
            _ZetaRadian(FLOAT_INIT), _KappaRadian(FLOAT_INIT),
            _EtaRadian(FLOAT_INIT) {};
        ZYZEulerAngles(Real zeta, Real kappa, Real eta) :
            _ZetaRadian(zeta), _KappaRadian(kappa), _EtaRadian(eta) {};
    };

    // angular step parameters structure
    struct ThetaAngles {
        Real _Theta1Degree;
        Real _Theta2Degree;
        Real _Theta3Degree;
        ThetaAngles() :
            _Theta1Degree(FLOAT_INIT), _Theta2Degree(FLOAT_INIT),
            _Theta3Degree(FLOAT_INIT) {};
        ThetaAngles(Real theta1, Real theta2, Real theta3) :
            _Theta1Degree(theta1), _Theta2Degree(theta2),
            _Theta3Degree(theta3) {};
    };

    // ZYZ Euler angles / angular step parameters conversion functions
    ThetaAngles ZYZEuler_2_Theta(const ZYZEulerAngles& euler_angles);
    ZYZEulerAngles Theta_2_ZYZEuler(const ThetaAngles& theta_angles);

    // step rotation matrix
    Matrix3 step_rotation_matrix(const ZYZEulerAngles& euler_angles);

    // step frame rotation matrix
    Matrix3 step_frame_rotation_matrix(const ZYZEulerAngles& euler_angles);

    // step Lambda vectors
    std::vector<Vector3> step_Lambda_vectors(const ZYZEulerAngles&
                                             euler_angles);
    Vector3 step_Lambda1_vector(const ZYZEulerAngles& euler_angles);
    Vector3 step_Lambda2_vector(const ZYZEulerAngles& euler_angles);
    Vector3 step_Lambda3_vector(const ZYZEulerAngles& euler_angles);

    // step Lambda vectors derivatives
    std::vector<Matrix3> step_dLambda_matrices(const ZYZEulerAngles&
                                               euler_angles);
    Matrix3 step_dLambda1_matrix(const ZYZEulerAngles& euler_angles);
    Matrix3 step_dLambda2_matrix(const ZYZEulerAngles& euler_angles);
    Matrix3 step_dLambda3_matrix(const ZYZEulerAngles& euler_angles);

    // step Xi matrix
    Matrix3 step_Xi_matrix(const ZYZEulerAngles& euler_angles);

    // step Omega matrix
    // this corresponds to the inverse of Xi
    Matrix3 step_Omega_matrix(const ZYZEulerAngles& euler_angles);

    // Xi matrix derivatives
    std::vector<Matrix3> step_dXi_matrices(const ZYZEulerAngles& euler_angles);
    Matrix3 step_Xid1_matrix(const ZYZEulerAngles& euler_angles);
    Matrix3 step_Xid2_matrix(const ZYZEulerAngles& euler_angles);
    Matrix3 step_Xid3_matrix(const ZYZEulerAngles& euler_angles);

    // Omega matrix derivatives
    Matrix3 step_Omegad1_matrix(const ZYZEulerAngles& euler_angles);
    Matrix3 step_Omegad2_matrix(const ZYZEulerAngles& euler_angles);
    Matrix3 step_Omegad3_matrix(const ZYZEulerAngles& euler_angles);


}


#endif  // emDNA_BpGeometryFunctions_h
