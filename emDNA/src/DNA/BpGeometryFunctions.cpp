// BpGeometryFunctions
// Nicolas Clauvelin


#include "DNA/BpStepDofs.h"
#include "DNA/BpStepParams.h"
#include "DNA/BpGeometryFunctions.h"


namespace BpGeometryFunctions {


    // ZYZ Euler angles / angular step parameters conversion
    ThetaAngles ZYZEuler_2_Theta(const ZYZEulerAngles& euler_angles) {

        // gamma angle
        const Real GammaRadian((euler_angles._EtaRadian-
                                euler_angles._ZetaRadian)/Real(2));

        // theta angles
        return ThetaAngles((euler_angles._KappaRadian*
                            std::sin(GammaRadian))*RAD_2_DEG,
                           (euler_angles._KappaRadian*
                            std::cos(GammaRadian))*RAD_2_DEG,
                           (euler_angles._ZetaRadian+
                            euler_angles._EtaRadian)*RAD_2_DEG);

    };
    ZYZEulerAngles Theta_2_ZYZEuler(const ThetaAngles& theta_angles) {

        ZYZEulerAngles euler_angles;

        // convert to radians
        const Real theta1Radians = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radians = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta3Radians = theta_angles._Theta3Degree*DEG_2_RAD;

        // Euler angles
        const Real GammaRadian = std::atan2(theta1Radians, theta2Radians);
        euler_angles._ZetaRadian = (theta3Radians/Real(2))-GammaRadian;
        euler_angles._KappaRadian = std::hypot(theta1Radians, theta2Radians);
        euler_angles._EtaRadian = (theta3Radians/Real(2))+GammaRadian;

        return euler_angles;

    };


    // step rotation matrix
    Matrix3 step_rotation_matrix(const ZYZEulerAngles& euler_angles) {

        // cosine and sine
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real cosKappa = std::cos(euler_angles._KappaRadian);
        const Real sinKappa = std::sin(euler_angles._KappaRadian);
        const Real cosEta = std::cos(euler_angles._EtaRadian);
        const Real sinEta = std::sin(euler_angles._EtaRadian);

        // array of values
        Real matelems[9];
        matelems[0] = cosEta*cosKappa*cosZeta - sinEta*sinZeta;
        matelems[1] = -(cosKappa*cosZeta*sinEta) - cosEta*sinZeta;
        matelems[2] = cosZeta*sinKappa;
        matelems[3] = cosZeta*sinEta + cosEta*cosKappa*sinZeta;
        matelems[4] = cosEta*cosZeta - cosKappa*sinEta*sinZeta;
        matelems[5] = sinKappa*sinZeta;
        matelems[6] = -(cosEta*sinKappa);
        matelems[7] = sinEta*sinKappa;
        matelems[8] = cosKappa;
        
        Matrix3 rotmat;
        rotmat << Array<Real>(9, matelems);

        return rotmat;

    };


    // step frame rotation matrix
    Matrix3 step_frame_rotation_matrix(const ZYZEulerAngles& euler_angles) {

        // gamma angle
        const Real GammaRadian =
            (euler_angles._EtaRadian-euler_angles._ZetaRadian)/Real(2);

        // cosines and sines
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2));
        const Real cosGamma = std::cos(GammaRadian);
        const Real sinGamma = std::sin(GammaRadian);

        // array of values
        Real matelems[9];
        matelems[0] = cosGamma*cosKappaHalf*cosZeta - sinGamma*sinZeta;
        matelems[1] = -(cosKappaHalf*cosZeta*sinGamma) - cosGamma*sinZeta;
        matelems[2] = cosZeta*sinKappaHalf;
        matelems[3] = cosZeta*sinGamma + cosGamma*cosKappaHalf*sinZeta;
        matelems[4] = cosGamma*cosZeta - cosKappaHalf*sinGamma*sinZeta;
        matelems[5] = sinKappaHalf*sinZeta;
        matelems[6] = -(cosGamma*sinKappaHalf);
        matelems[7] = sinGamma*sinKappaHalf;
        matelems[8] = cosKappaHalf;

        Matrix3 rotmat;
        rotmat << Array<Real>(9, matelems);

        return rotmat;

    };


    // step Lambda vectors
    std::vector<Vector3> step_Lambda_vectors(const ZYZEulerAngles&
                                             euler_angles) {

        // precomputed quantities
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        // vectors
        Vector3 Lambda1, Lambda2, Lambda3;

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            Lambda1[X] = -Real(1)/Real(2);  // Lambda1
            Lambda1[Y] = FLOAT_INIT;
            Lambda1[Z] = FLOAT_INIT;
            Lambda2[X] = FLOAT_INIT;         // Lambda2
            Lambda2[Y] = -Real(1)/Real(2);
            Lambda2[Z] = FLOAT_INIT;
            Lambda3[X] = FLOAT_INIT;         // Lambda3
            Lambda3[Y] = FLOAT_INIT;
            Lambda3[Z] = -Real(1)/Real(2);
        }

        // regular case
        else {
            // Lambda1
            Lambda1[X] = (-theta1RadianSq/(Real(2)*kappaRadianSq) -
                          (sinKappaHalf*theta2RadianSq)/
                          (kappaRadianSq*kappaRadian));
            Lambda1[Y] = (-((kappaRadian - Real(2)*sinKappaHalf)*
                            theta1Radian*theta2Radian)/
                          (Real(2)*kappaRadianSq*kappaRadian));
            Lambda1[Z] = (((-1 + cosKappaHalf)*theta2Radian)/kappaRadianSq);
            // Lambda2
            Lambda2[X] = (-((kappaRadian - Real(2)*sinKappaHalf)*
                            theta1Radian*theta2Radian)/
                          (Real(2)*kappaRadianSq*kappaRadian));
            Lambda2[Y] = (-(sinKappaHalf*theta1RadianSq)/
                          (kappaRadianSq*kappaRadian)
                          - theta2RadianSq/(Real(2)*kappaRadianSq));
            Lambda2[Z] = (-(((-1 + cosKappaHalf)*theta1Radian)/kappaRadianSq));
            // Lambda3
            Lambda3[X] = (sinKappaHalf*theta2Radian)/(Real(2)*kappaRadian);
            Lambda3[Y] = -(sinKappaHalf*theta1Radian)/(Real(2)*kappaRadian);
            Lambda3[Z] = -cosKappaHalf/Real(2);
        };

        std::vector<Vector3> Lambda_vectors;
        Lambda_vectors.reserve(3);
        Lambda_vectors.push_back(Lambda1);
        Lambda_vectors.push_back(Lambda2);
        Lambda_vectors.push_back(Lambda3);

        return Lambda_vectors;

    };
    Vector3 step_Lambda1_vector(const ZYZEulerAngles& euler_angles) {

        // cosine and sine
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        // vector
        Vector3 Lambda1;

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            Lambda1[X] = -Real(1)/Real(2);
            Lambda1[Y] = FLOAT_INIT;
            Lambda1[Z] = FLOAT_INIT;
        }

        // regular case
        else {
            Lambda1[X] = (-theta1RadianSq/(Real(2)*kappaRadianSq) -
                          (sinKappaHalf*theta2RadianSq)/
                          (kappaRadianSq*kappaRadian));
            Lambda1[Y] = (-((kappaRadian - Real(2)*sinKappaHalf)*
                            theta1Radian*theta2Radian)/
                          (Real(2)*kappaRadianSq*kappaRadian));
            Lambda1[Z] = (((-1 + cosKappaHalf)*theta2Radian)/kappaRadianSq);
        };

        return Lambda1;

    };
    Vector3 step_Lambda2_vector(const ZYZEulerAngles& euler_angles) {

        // cosine and sine
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        // vector
        Vector3 Lambda2;

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            Lambda2[X] = FLOAT_INIT;
            Lambda2[Y] = -Real(1)/Real(2);
            Lambda2[Z] = FLOAT_INIT;
        }

        // regular case
        else {
            Lambda2[X] = (-((kappaRadian - Real(2)*sinKappaHalf)*
                            theta1Radian*theta2Radian)/
                          (Real(2)*kappaRadianSq*kappaRadian));
            Lambda2[Y] = (-(sinKappaHalf*theta1RadianSq)/
                            (kappaRadianSq*kappaRadian)
                          - theta2RadianSq/(Real(2)*kappaRadianSq));
            Lambda2[Z] = (-(((-1 + cosKappaHalf)*theta1Radian)/kappaRadianSq));
        };
        
        return Lambda2;

    };
    Vector3 step_Lambda3_vector(const ZYZEulerAngles& euler_angles) {

        // cosine and sine
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2));
        const Real kappaRadian = euler_angles._KappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;

        // vector
        Vector3 Lambda3;

        // check for degenerated case
        if (euler_angles._KappaRadian*euler_angles._KappaRadian < ZERO_EPS) {
            Lambda3[X] = FLOAT_INIT;
            Lambda3[Y] = FLOAT_INIT;
            Lambda3[Z] = -Real(1)/Real(2);
        }

        // regular case
        else {
            Lambda3[X] = (sinKappaHalf*theta2Radian)/(Real(2)*kappaRadian);
            Lambda3[Y] = -(sinKappaHalf*theta1Radian)/(Real(2)*kappaRadian);
            Lambda3[Z] = -cosKappaHalf/Real(2);
        };
        
        return Lambda3;
        
    };


    // step Lambda vectors derivatives
    std::vector<Matrix3> step_dLambda_matrices(const ZYZEulerAngles&
                                               euler_angles) {

        // precomputed quantities
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2.0));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2.0));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        Matrix3 dLambda1, dLambda2, dLambda3;

        // check for degenerated case
        if (euler_angles._KappaRadian*euler_angles._KappaRadian < ZERO_EPS) {
            dLambda1(2,1) = -0.125;
            dLambda2(2,0) = 0.125;
            dLambda3(0,1) = 0.25;
            dLambda3(1,0) = -0.25;
        }

        // regular case
        else {

            // dLambda1
            dLambda1(0,0) =
            -theta1Radian*theta2RadianSq/(kappaRadianSq*kappaRadianSq)*
            (Real(1.0)+Real(0.5)*cosKappaHalf
             -Real(3.0)*sinKappaHalf/kappaRadian);
            dLambda1(0,1) =
            theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (theta1RadianSq-Real(0.5)*theta2RadianSq*cosKappaHalf
             -Real(2.0)*(theta1RadianSq-Real(0.5)*theta2RadianSq)*
             sinKappaHalf/kappaRadian);
            dLambda1(1,0) =
            theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta1RadianSq*cosKappaHalf+
             (theta2RadianSq-Real(2.0)*theta1RadianSq)*sinKappaHalf/
             kappaRadian);
            dLambda1(1,1) =
            theta1Radian/(kappaRadianSq*kappaRadianSq)*
            (-Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta2RadianSq*cosKappaHalf+
             (theta1RadianSq-Real(2.0)*theta2RadianSq)*sinKappaHalf/
             kappaRadian);
            dLambda1(2,0) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (-Real(2.0)*(cosKappaHalf-Real(1.0))/kappaRadianSq
             -Real(0.5)*sinKappaHalf/kappaRadian);
            dLambda1(2,1) =
            Real(1.0)/kappaRadianSq*
            ((theta1RadianSq-theta2RadianSq)*(cosKappaHalf-Real(1.0))/
             kappaRadianSq
             -Real(0.5)*theta2RadianSq*sinKappaHalf/kappaRadian);

            // dLambda2
            dLambda2(0,0) =
            theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta1RadianSq*cosKappaHalf+
             (theta2RadianSq-Real(2.0)*theta1RadianSq)*sinKappaHalf/
             kappaRadian);
            dLambda2(0,1) =
            theta1Radian/(kappaRadianSq*kappaRadianSq)*
            (-Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta2RadianSq*cosKappaHalf+
             (theta1RadianSq-Real(2.0)*theta2RadianSq)*sinKappaHalf/
             kappaRadian);
            dLambda2(1,0) =
            theta1Radian/(kappaRadianSq*kappaRadianSq)*
            (theta2RadianSq-Real(0.5)*theta1RadianSq*cosKappaHalf
             +(theta1RadianSq-Real(2.0)*theta2RadianSq)*sinKappaHalf/
             kappaRadian);
            dLambda2(1,1) =
            theta1RadianSq*theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (-Real(1.0)-Real(0.5)*cosKappaHalf+Real(3.0)*sinKappaHalf/
             kappaRadian);
            dLambda2(2,0) =
            Real(1.0)/kappaRadianSq*
            ((theta1RadianSq-theta2RadianSq)*(cosKappaHalf-Real(1.0))/
             kappaRadianSq
             +Real(0.5)*theta1RadianSq*sinKappaHalf/kappaRadian);
            dLambda2(2,1) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (Real(2.0)*(cosKappaHalf-Real(1.0))/kappaRadianSq
             +Real(0.5)*sinKappaHalf/kappaRadian);

            // dLambda3
            dLambda3(0,0) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (Real(0.25)*cosKappaHalf-Real(0.5)*sinKappaHalf/kappaRadian);
            dLambda3(0,1) =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (Real(0.25)*kappaRadian*
             theta2RadianSq*(theta1RadianSq+theta2RadianSq)*cosKappaHalf
             +Real(0.5)*theta1RadianSq*(theta1RadianSq+theta2RadianSq)
             *sinKappaHalf);
            dLambda3(1,0) =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (-Real(0.25)*kappaRadian*
             theta1RadianSq*(theta1RadianSq+theta2RadianSq)*cosKappaHalf
             -Real(0.5)*theta2RadianSq*(theta1RadianSq+theta2RadianSq)
             *sinKappaHalf);
            dLambda3(1,1) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (-Real(0.25)*cosKappaHalf+Real(0.5)*sinKappaHalf/kappaRadian);
            dLambda3(2,0) =
            theta1Radian*sinKappaHalf/(Real(4.0)*kappaRadian);;
            dLambda3(2,1) =
            theta2Radian*sinKappaHalf/(Real(4.0)*kappaRadian);

        };

        std::vector<Matrix3> dLamda_mats;
        dLamda_mats.push_back(dLambda1);
        dLamda_mats.push_back(dLambda2);
        dLamda_mats.push_back(dLambda3);

        return dLamda_mats;

    };
    Matrix3 step_dLambda1_matrix(const ZYZEulerAngles& euler_angles) {

        // precomputed quantities
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2.0));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2.0));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        Matrix3 dLambda1;

        // check for degenerated case
        if (euler_angles._KappaRadian*euler_angles._KappaRadian < ZERO_EPS) {
            dLambda1(2,1) = -0.125;
        }

        // regular case
        else {

            // entry 0,0
            dLambda1(0,0) =
            -theta1Radian*theta2RadianSq/(kappaRadianSq*kappaRadianSq)*
            (Real(1.0)+Real(0.5)*cosKappaHalf
             -Real(3.0)*sinKappaHalf/kappaRadian);

            // entry 0,1
            dLambda1(0,1) =
            theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (theta1RadianSq-Real(0.5)*theta2RadianSq*cosKappaHalf
             -Real(2.0)*(theta1RadianSq-Real(0.5)*theta2RadianSq)*
             sinKappaHalf/kappaRadian);

            // entry 1,0
            dLambda1(1,0) =
            theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta1RadianSq*cosKappaHalf+
             (theta2RadianSq-Real(2.0)*theta1RadianSq)*sinKappaHalf/
             kappaRadian);

            // entry 1,1
            dLambda1(1,1) =
            theta1Radian/(kappaRadianSq*kappaRadianSq)*
            (-Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta2RadianSq*cosKappaHalf+
             (theta1RadianSq-Real(2.0)*theta2RadianSq)*sinKappaHalf/
             kappaRadian);

            // entry 2,0
            dLambda1(2,0) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (-Real(2.0)*(cosKappaHalf-Real(1.0))/kappaRadianSq
             -Real(0.5)*sinKappaHalf/kappaRadian);

            // entry 2,1
            dLambda1(2,1) =
            Real(1.0)/kappaRadianSq*
            ((theta1RadianSq-theta2RadianSq)*(cosKappaHalf-Real(1.0))/
             kappaRadianSq
             -Real(0.5)*theta2RadianSq*sinKappaHalf/kappaRadian);

        };

        return dLambda1;

    };
    Matrix3 step_dLambda2_matrix(const ZYZEulerAngles& euler_angles) {

        // precomputed quantities
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2.0));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2.0));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        Matrix3 dLambda2;

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            dLambda2(2,0) = 0.125;
        }

        // regular case
        else {

            // entry 0,0
            dLambda2(0,0) =
            theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta1RadianSq*cosKappaHalf+
             (theta2RadianSq-Real(2.0)*theta1RadianSq)*sinKappaHalf/
             kappaRadian);

            // entry 0,1
            dLambda2(0,1) =
            theta1Radian/(kappaRadianSq*kappaRadianSq)*
            (-Real(0.5)*(theta1RadianSq-theta2RadianSq)+
             Real(0.5)*theta2RadianSq*cosKappaHalf+
             (theta1RadianSq-Real(2.0)*theta2RadianSq)*sinKappaHalf/
             kappaRadian);

            // entry 1,0
            dLambda2(1,0) =
            theta1Radian/(kappaRadianSq*kappaRadianSq)*
            (theta2RadianSq-Real(0.5)*theta1RadianSq*cosKappaHalf
             +(theta1RadianSq-Real(2.0)*theta2RadianSq)*sinKappaHalf/
             kappaRadian);

            // entry 1,1
            dLambda2(1,1) =
            theta1RadianSq*theta2Radian/(kappaRadianSq*kappaRadianSq)*
            (-Real(1.0)-Real(0.5)*cosKappaHalf+Real(3.0)*sinKappaHalf/
             kappaRadian);

            // entry 2,0
            dLambda2(2,0) =
            Real(1.0)/kappaRadianSq*
            ((theta1RadianSq-theta2RadianSq)*(cosKappaHalf-Real(1.0))/
             kappaRadianSq
             +Real(0.5)*theta1RadianSq*sinKappaHalf/kappaRadian);

            // entry 2,1
            dLambda2(2,1) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (Real(2.0)*(cosKappaHalf-Real(1.0))/kappaRadianSq
             +Real(0.5)*sinKappaHalf/kappaRadian);
            
        };
        
        return dLambda2;
        
    };
    Matrix3 step_dLambda3_matrix(const ZYZEulerAngles& euler_angles) {

        // precomputed quantities
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2.0));
        const Real sinKappaHalf = std::sin(euler_angles._KappaRadian/Real(2.0));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;

        Matrix3 dLambda3;

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            dLambda3(0,1) = 0.25;
            dLambda3(1,0) = -0.25;
        }

        // regular case
        else {

            // entry 0,0
            dLambda3(0,0) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (Real(0.25)*cosKappaHalf-Real(0.5)*sinKappaHalf/kappaRadian);

            // entry 0,1
            dLambda3(0,1) =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (Real(0.25)*kappaRadian*
             theta2RadianSq*(theta1RadianSq+theta2RadianSq)*cosKappaHalf
             +Real(0.5)*theta1RadianSq*(theta1RadianSq+theta2RadianSq)
             *sinKappaHalf);

            // entry 1,0
            dLambda3(1,0) =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (-Real(0.25)*kappaRadian*
            theta1RadianSq*(theta1RadianSq+theta2RadianSq)*cosKappaHalf
             -Real(0.5)*theta2RadianSq*(theta1RadianSq+theta2RadianSq)
             *sinKappaHalf);

            // entry 1,1
            dLambda3(1,1) =
            theta1Radian*theta2Radian/kappaRadianSq*
            (-Real(0.25)*cosKappaHalf+Real(0.5)*sinKappaHalf/kappaRadian);

            // entry 2,0
            dLambda3(2,0) =
            theta1Radian*sinKappaHalf/(Real(4.0)*kappaRadian);;

            // entry 2,1
            dLambda3(2,1) =
            theta2Radian*sinKappaHalf/(Real(4.0)*kappaRadian);
            
        };
        
        return dLambda3;

    };
    

    // step Xi matrix
    Matrix3 step_Xi_matrix(const ZYZEulerAngles& euler_angles) {

        // cosines and sines
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real cosKappa = std::cos(euler_angles._KappaRadian);
        const Real sinKappa = std::sin(euler_angles._KappaRadian);
        const Real cosKappaHalf = std::cos(euler_angles._KappaRadian/Real(2));
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;

        // array of values
        Real matelems[9];

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            matelems[0] = cosZeta;
            matelems[1] = -sinZeta;
            matelems[2] = FLOAT_INIT;
            matelems[3] = sinZeta;
            matelems[4] = cosZeta;
            matelems[5] = FLOAT_INIT;
            matelems[6] = FLOAT_INIT;
            matelems[7] = FLOAT_INIT;
            matelems[8] = Real(1);
        }

        // regular case
        else {
            matelems[0] = (-kappaRadian*sinZeta*theta1Radian +
                           cosZeta*sinKappa*theta2Radian)/kappaRadianSq;
            matelems[1] = -(cosZeta*sinKappa*theta1Radian +
                            kappaRadian*sinZeta*theta2Radian)/kappaRadianSq;
            matelems[2] = (cosZeta*sinKappa)/Real(2);
            matelems[3] = (cosZeta*kappaRadian*theta1Radian +
                           sinKappa*sinZeta*theta2Radian)/kappaRadianSq;
            matelems[4] = (-sinKappa*sinZeta*theta1Radian +
                           cosZeta*kappaRadian*theta2Radian)/kappaRadianSq;
            matelems[5] = (sinKappa*sinZeta)/Real(2);
            matelems[6] = ((-1 + cosKappa)*theta2Radian)/kappaRadianSq;
            matelems[7] = (theta1Radian - cosKappa*theta1Radian)/kappaRadianSq;
            matelems[8] = cosKappaHalf*cosKappaHalf;
        };

        Matrix3 ximat;
        ximat << Array<Real>(9, matelems);

        return ximat;

    };


    // step Sigma matrix
    Matrix3 step_Omega_matrix(const ZYZEulerAngles& euler_angles) {

        // cosine and sine
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real tanKappaHalf = std::tan(euler_angles._KappaRadian/Real(2));
        const Real kappaRadian = euler_angles._KappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;

        // array of values
        Real matelems[9];

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {
            matelems[0] = cosZeta;
            matelems[1] = -sinZeta;
            matelems[2] = FLOAT_INIT;
            matelems[3] = sinZeta;
            matelems[4] = cosZeta;
            matelems[5] = FLOAT_INIT;
            matelems[6] = FLOAT_INIT;
            matelems[7] = FLOAT_INIT;
            matelems[8] = Real(1);
        }

        // regular case
        else {
            matelems[0] = (-((sinZeta*theta1Radian)/kappaRadian)
                           + (cosZeta*theta2Radian)/(Real(2)*tanKappaHalf));
            matelems[1] = (-(cosZeta*theta1Radian)/(Real(2)*tanKappaHalf)
                           - (sinZeta*theta2Radian)/kappaRadian);
            matelems[2] = cosZeta*tanKappaHalf;
            matelems[3] = ((cosZeta*theta1Radian)/kappaRadian +
                           (sinZeta*theta2Radian)/(Real(2)*tanKappaHalf));
            matelems[4] = (-(sinZeta*theta1Radian)/(Real(2)*tanKappaHalf) +
                           (cosZeta*theta2Radian)/kappaRadian);
            matelems[5] = sinZeta*tanKappaHalf;
            matelems[6] = -theta2Radian/Real(2);
            matelems[7] = theta1Radian/Real(2);
            matelems[8] = Real(1);
        };

        // we perform a transpose cause we compute the inverse transpose
        // at this point
        Matrix3 omegamat;
        omegamat << Array<Real>(9, matelems);
        omegamat.transpose();
        
        return omegamat;

    };


    // Xi matrix derivatives
    std::vector<Matrix3> step_dXi_matrices(const ZYZEulerAngles& euler_angles) {

        std::vector<Matrix3> dXi_mats;
        dXi_mats.push_back(step_Xid1_matrix(euler_angles));
        dXi_mats.push_back(step_Xid2_matrix(euler_angles));
        dXi_mats.push_back(step_Xid3_matrix(euler_angles));

        return dXi_mats;

    };
    Matrix3 step_Xid1_matrix(const ZYZEulerAngles& euler_angles) {

        // precomputed quantities
        const Real cosKappa = std::cos(euler_angles._KappaRadian);
        const Real sinKappa = std::sin(euler_angles._KappaRadian);
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta3Radian = theta_angles._Theta3Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;
        const Real cosTheta3Half = std::cos(Real(0.5)*theta3Radian);
        const Real sinTheta3Half = std::sin(Real(0.5)*theta3Radian);

        Real dXi_1_elems[9];

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {

            dXi_1_elems[0] = FLOAT_INIT;
            dXi_1_elems[1] = FLOAT_INIT;
            dXi_1_elems[2] = Real(0.5)*sinTheta3Half;
            dXi_1_elems[3] = FLOAT_INIT;
            dXi_1_elems[4] = FLOAT_INIT;
            dXi_1_elems[5] = -Real(0.5)*cosTheta3Half;
            dXi_1_elems[6] = FLOAT_INIT;
            dXi_1_elems[7] = Real(0.5);
            dXi_1_elems[8] = FLOAT_INIT;

        }

        // regular case
        else {

            dXi_1_elems[0] =
            theta2Radian/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (theta1Radian*kappaRadianSq*cosZeta*(Real(1.0)+cosKappa)-
             theta2Radian*kappaRadianSq*sinZeta+
             (-Real(3.0)*theta1Radian*theta2Radian*cosTheta3Half+
              (-Real(2.0)*theta1RadianSq+theta2RadianSq)*sinTheta3Half)*sinKappa);

            dXi_1_elems[1] =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (kappaRadianSq*cosZeta*(theta2RadianSq-theta1RadianSq*cosKappa)+
             theta1Radian*theta2Radian*kappaRadianSq*sinZeta+
             (theta2Radian*(-theta2RadianSq+Real(2.0)*theta1RadianSq)*cosTheta3Half+
              theta1Radian*(theta1RadianSq-Real(2.0)*theta2RadianSq)*sinTheta3Half)*
             sinKappa);

            dXi_1_elems[2] =
            Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (theta1Radian*kappaRadian*cosZeta*cosKappa+
             theta2Radian*sinZeta*sinKappa);

            dXi_1_elems[3] =
            theta2Radian/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (theta2Radian*kappaRadianSq*cosZeta+
             theta1Radian*kappaRadianSq*(Real(1.0)+cosKappa)*sinZeta+
             ((Real(2.0)*theta1RadianSq-theta2RadianSq)*cosTheta3Half+
              -Real(3.0)*theta1Radian*theta2Radian*sinTheta3Half)*sinKappa);

            dXi_1_elems[4] =
            -Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (-kappaRadianSq*(theta2RadianSq-theta1RadianSq*cosKappa)*sinZeta+
             theta1Radian*theta2Radian*kappaRadianSq*cosZeta+
             (theta1Radian*(theta1RadianSq-Real(2.0)*theta2RadianSq)*cosTheta3Half+
              theta2Radian*(theta2RadianSq-Real(2.0)*theta1RadianSq)*sinTheta3Half)*
             sinKappa);

            dXi_1_elems[5] =
            Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (theta1Radian*kappaRadian*cosKappa*sinZeta-
             theta2Radian*cosZeta*sinKappa);

            dXi_1_elems[6] =
            -(theta1Radian*theta2Radian)/(kappaRadianSq*kappaRadianSq)*
            (-Real(2.0)+Real(2.0)*cosKappa+kappaRadian*sinKappa);

            dXi_1_elems[7] =
            Real(1.0)/(kappaRadianSq*kappaRadianSq)*
            ((theta2RadianSq-theta1RadianSq)+
             (theta1RadianSq-theta2RadianSq)*cosKappa+
             theta1RadianSq*kappaRadian*sinKappa);

            dXi_1_elems[8] = -theta1Radian*sinKappa/(Real(2.0)*kappaRadian);

        };

        Matrix3 dXimat;
        dXimat << Array<Real>(9, dXi_1_elems);

        return dXimat;

    };
    Matrix3 step_Xid2_matrix(const ZYZEulerAngles& euler_angles) {

        // precomputed quantities
        const Real cosKappa = std::cos(euler_angles._KappaRadian);
        const Real sinKappa = std::sin(euler_angles._KappaRadian);
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta3Radian = theta_angles._Theta3Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;
        const Real cosTheta3Half = std::cos(Real(0.5)*theta3Radian);
        const Real sinTheta3Half = std::sin(Real(0.5)*theta3Radian);

        Real dXi_2_elems[9];

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {

            dXi_2_elems[0] = FLOAT_INIT;
            dXi_2_elems[1] = FLOAT_INIT;
            dXi_2_elems[2] = Real(0.5)*cosTheta3Half;
            dXi_2_elems[3] = FLOAT_INIT;
            dXi_2_elems[4] = FLOAT_INIT;
            dXi_2_elems[5] = Real(0.5)*sinTheta3Half;
            dXi_2_elems[6] = -Real(0.5);
            dXi_2_elems[7] = FLOAT_INIT;
            dXi_2_elems[8] = FLOAT_INIT;

        }

        // regular case
        else {

            dXi_2_elems[0] =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (-kappaRadianSq*cosZeta*(theta1RadianSq-theta2RadianSq*cosKappa)+
             theta1Radian*theta2Radian*kappaRadianSq*sinZeta+
             (theta2Radian*(-theta2RadianSq+Real(2.0)*theta1RadianSq)*cosTheta3Half+
              theta1Radian*(theta1RadianSq-Real(2.0)*theta2RadianSq)*sinTheta3Half)*
             sinKappa);

            dXi_2_elems[1] =
            theta1Radian/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (-theta2Radian*kappaRadianSq*cosZeta*(Real(1.0)+cosKappa)-
             theta1Radian*kappaRadianSq*sinZeta+
             (Real(3.0)*theta1Radian*theta2Radian*sinTheta3Half+
              (Real(2.0)*theta2RadianSq-theta1RadianSq)*cosTheta3Half)*sinKappa);

            dXi_2_elems[2] =
            Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (theta2Radian*kappaRadian*cosZeta*cosKappa-
             theta1Radian*sinZeta*sinKappa);

            dXi_2_elems[3] =
            Real(1.0)/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (-kappaRadianSq*(theta1RadianSq-theta2RadianSq*cosKappa)*sinZeta
             -theta1Radian*theta2Radian*kappaRadianSq*cosZeta-
             (theta1Radian*(theta1RadianSq-Real(2.0)*theta2RadianSq)*cosTheta3Half+
              theta2Radian*(theta2RadianSq-Real(2.0)*theta1RadianSq)*sinTheta3Half)*
             sinKappa);

            dXi_2_elems[4] =
            theta1Radian/(kappaRadianSq*kappaRadianSq*kappaRadian)*
            (theta1Radian*kappaRadianSq*cosZeta
             -theta2Radian*kappaRadianSq*(Real(1.0)+cosKappa)*sinZeta-
             ((-Real(2.0)*theta2RadianSq+theta1RadianSq)*sinTheta3Half+
              Real(3.0)*theta1Radian*theta2Radian*cosTheta3Half)*sinKappa);

            dXi_2_elems[5] =
            Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (theta2Radian*kappaRadian*cosKappa*sinZeta+
             theta1Radian*cosZeta*sinKappa);

            dXi_2_elems[6] =
            Real(1.0)/(kappaRadianSq*kappaRadianSq)*
            ((theta2RadianSq-theta1RadianSq)+
             (theta1RadianSq-theta2RadianSq)*cosKappa-
             theta2RadianSq*kappaRadian*sinKappa);

            dXi_2_elems[7] =
            (theta1Radian*theta2Radian)/(kappaRadianSq*kappaRadianSq)*
            (-Real(2.0)+Real(2.0)*cosKappa+kappaRadian*sinKappa);

            dXi_2_elems[8] = -theta2Radian*sinKappa/(Real(2.0)*kappaRadian);
            
        };
        
        Matrix3 dXimat;
        dXimat << Array<Real>(9, dXi_2_elems);
        
        return dXimat;

    };
    Matrix3 step_Xid3_matrix(const ZYZEulerAngles& euler_angles) {

        // precomputed quantities
        const Real sinKappa = std::sin(euler_angles._KappaRadian);
        const Real cosZeta = std::cos(euler_angles._ZetaRadian);
        const Real sinZeta = std::sin(euler_angles._ZetaRadian);
        const Real kappaRadian = euler_angles._KappaRadian;
        const Real kappaRadianSq = kappaRadian*kappaRadian;

        // theta angles in radians
        const ThetaAngles theta_angles = ZYZEuler_2_Theta(euler_angles);
        const Real theta1Radian = theta_angles._Theta1Degree*DEG_2_RAD;
        const Real theta2Radian = theta_angles._Theta2Degree*DEG_2_RAD;
        const Real theta3Radian = theta_angles._Theta3Degree*DEG_2_RAD;
        const Real theta1RadianSq = theta1Radian*theta1Radian;
        const Real theta2RadianSq = theta2Radian*theta2Radian;
        const Real cosTheta3Half = std::cos(Real(0.5)*theta3Radian);
        const Real sinTheta3Half = std::sin(Real(0.5)*theta3Radian);

        Real dXi_3_elems[9];

        // check for degenerated case
        if (euler_angles._KappaRadian < ZERO_EPS) {

            dXi_3_elems[0] = -Real(0.5)*sinTheta3Half;
            dXi_3_elems[1] = -Real(0.5)*cosTheta3Half;
            dXi_3_elems[2] = FLOAT_INIT;
            dXi_3_elems[3] = Real(0.5)*sinTheta3Half;
            dXi_3_elems[4] = -Real(0.5)*sinTheta3Half;
            dXi_3_elems[5] = FLOAT_INIT;
            dXi_3_elems[6] = FLOAT_INIT;
            dXi_3_elems[7] = FLOAT_INIT;
            dXi_3_elems[8] = FLOAT_INIT;

        }

        // regular case
        else {

            dXi_3_elems[0] =
            -Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (theta1Radian*theta2Radian*cosTheta3Half+
             theta1RadianSq*sinTheta3Half+
             theta2Radian*sinZeta*sinKappa);

            dXi_3_elems[1] =
            Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (-theta2RadianSq*cosTheta3Half
             -theta1Radian*theta2Radian*sinTheta3Half
             +theta1Radian*sinZeta*sinKappa);

            dXi_3_elems[2] = -Real(0.25)*sinZeta*sinKappa;

            dXi_3_elems[3] =
            Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (theta1RadianSq*cosTheta3Half
             -theta1Radian*theta2Radian*sinTheta3Half
             +theta2Radian*cosZeta*sinKappa);

            dXi_3_elems[4] =
            -Real(1.0)/(Real(2.0)*kappaRadianSq)*
            (-theta1Radian*theta2Radian*cosTheta3Half+
             theta2RadianSq*sinTheta3Half+
             theta1Radian*cosZeta*sinKappa);

            dXi_3_elems[5] = Real(0.25)*cosZeta*sinKappa;

            dXi_3_elems[6] = FLOAT_INIT;
            dXi_3_elems[7] = FLOAT_INIT;
            dXi_3_elems[8] = FLOAT_INIT;
            
        };
        
        Matrix3 dXimat;
        dXimat << Array<Real>(9, dXi_3_elems);
        
        return dXimat;

    };


}

