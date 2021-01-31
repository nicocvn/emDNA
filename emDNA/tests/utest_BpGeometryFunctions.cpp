// BpGeometricFunctions unit tests
// Nicolas Clauvelin


#include <utest_BpGeometryFunctions.h>


namespace {

    
    // ZYZEulerAngles <-> ThetaAngles
	TEST_F(BpGeometryFunctionsTest, ZYZEulerAnglesAndThetaAngles) {

        ThetaAngles th_angles = ZYZEuler_2_Theta(test_euler_angles);
        ZYZEulerAngles eu_angles = Theta_2_ZYZEuler(test_theta_angles);

        EXPECT_NEAR(test_theta_angles._Theta1Degree, th_angles._Theta1Degree,
                    ZERO_TOL);
        EXPECT_NEAR(test_theta_angles._Theta2Degree, th_angles._Theta2Degree,
                    ZERO_TOL);
        EXPECT_NEAR(test_theta_angles._Theta3Degree, th_angles._Theta3Degree,
                    ZERO_TOL);
        EXPECT_NEAR(test_euler_angles._ZetaRadian, eu_angles._ZetaRadian,
                    ZERO_TOL);
        EXPECT_NEAR(test_euler_angles._KappaRadian, eu_angles._KappaRadian,
                    ZERO_TOL);
        EXPECT_NEAR(test_euler_angles._EtaRadian, eu_angles._EtaRadian,
                    ZERO_TOL);

	};


    // Xi matrix
    TEST_F(BpGeometryFunctionsTest, XiMatrix) {

        Matrix3 ximat = step_Xi_matrix(test_euler_angles);
        Matrix3 dximat = ximat-mm_ximat;

        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                EXPECT_NEAR(mm_ximat(i,j), ximat(i,j), ZERO_EPS);

    };
    TEST_F(BpGeometryFunctionsTest, XiMatrixLimit) {

        ZYZEulerAngles angles = Theta_2_ZYZEuler(test_theta_angles_degen);
        Matrix3 ximat = step_Xi_matrix(angles);

        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                EXPECT_NEAR(mm_ximat_limit(i,j), ximat(i,j), ZERO_EPS);

    };


    // Xi matrix derivatives
    TEST_F(BpGeometryFunctionsTest, dXiMatrix) {

        // dXi matrices
        std::vector<Matrix3> dXi_mats =
        step_dXi_matrices(test_euler_angles);

        // checking
        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j) {
                EXPECT_NEAR(mm_dXi1(i,j), dXi_mats[0](i,j), ZERO_EPS);
                EXPECT_NEAR(mm_dXi2(i,j), dXi_mats[1](i,j), ZERO_EPS);
                EXPECT_NEAR(mm_dXi3(i,j), dXi_mats[2](i,j), ZERO_EPS);
            };

    };


    // Lambda vectors
    TEST_F(BpGeometryFunctionsTest, LambdaVectors) {

        Vector3 Lambda1 = step_Lambda1_vector(test_euler_angles);
        Vector3 Lambda2 = step_Lambda2_vector(test_euler_angles);
        Vector3 Lambda3 = step_Lambda3_vector(test_euler_angles);

        EXPECT_NEAR(mm_Lambda1[X], Lambda1[X], ZERO_EPS);
        EXPECT_NEAR(mm_Lambda1[Y], Lambda1[Y], ZERO_EPS);
        EXPECT_NEAR(mm_Lambda1[Z], Lambda1[Z], ZERO_EPS);

        EXPECT_NEAR(mm_Lambda2[X], Lambda2[X], ZERO_EPS);
        EXPECT_NEAR(mm_Lambda2[Y], Lambda2[Y], ZERO_EPS);
        EXPECT_NEAR(mm_Lambda2[Z], Lambda2[Z], ZERO_EPS);

        EXPECT_NEAR(mm_Lambda3[X], Lambda3[X], ZERO_EPS);
        EXPECT_NEAR(mm_Lambda3[Y], Lambda3[Y], ZERO_EPS);
        EXPECT_NEAR(mm_Lambda3[Z], Lambda3[Z], ZERO_EPS);

    };
    TEST_F(BpGeometryFunctionsTest, LambdaVectorsLimit) {

        ZYZEulerAngles angles = Theta_2_ZYZEuler(test_theta_angles_degen);

        Vector3 Lambda1 = step_Lambda1_vector(angles);
        Vector3 Lambda2 = step_Lambda2_vector(angles);
        Vector3 Lambda3 = step_Lambda3_vector(angles);

        EXPECT_NEAR(Real(-.5), Lambda1[X], ZERO_EPS);
        EXPECT_NEAR(0, Lambda1[Y], ZERO_EPS);
        EXPECT_NEAR(0, Lambda1[Z], ZERO_EPS);

        EXPECT_NEAR(0, Lambda2[X], ZERO_EPS);
        EXPECT_NEAR(Real(-.5), Lambda2[Y], ZERO_EPS);
        EXPECT_NEAR(0, Lambda2[Z], ZERO_EPS);

        EXPECT_NEAR(0, Lambda3[X], ZERO_EPS);
        EXPECT_NEAR(0, Lambda3[Y], ZERO_EPS);
        EXPECT_NEAR(Real(-.5), Lambda3[Z], ZERO_EPS);

    };


    // dLambda matrices
    TEST_F(BpGeometryFunctionsTest, dLambdaMatrices) {

        // dLambda matrices
        std::vector<Matrix3> dLambda_mats =
        step_dLambda_matrices(test_euler_angles);

        // checking
        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j) {
                EXPECT_NEAR(mm_dL1(i,j), dLambda_mats[0](i,j), ZERO_TOL);
                EXPECT_NEAR(mm_dL2(i,j), dLambda_mats[1](i,j), ZERO_TOL);
                EXPECT_NEAR(mm_dL3(i,j), dLambda_mats[2](i,j), ZERO_TOL);
            };

    };

    
    // Omega matrix
    // this matrix is such that its transpose is the inverse of Xi
    TEST_F(BpGeometryFunctionsTest, OmegaMatrix) {

        Matrix3 sigmamat = step_Omega_matrix(test_euler_angles);
        Matrix3 ximat = step_Xi_matrix(test_euler_angles);

        sigmamat *= ximat;
        sigmamat -= Matrix3::identity_matrix();

        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                EXPECT_NEAR(FLOAT_INIT, sigmamat(i,j), ZERO_EPS);
        
    };
    TEST_F(BpGeometryFunctionsTest, OmegasMatrixLimit) {

        ZYZEulerAngles angles = Theta_2_ZYZEuler(test_theta_angles_degen);
        Matrix3 sigmamat = step_Omega_matrix(angles);
        Matrix3 ximat = step_Xi_matrix(angles);

        sigmamat *= ximat;
        sigmamat -= Matrix3::identity_matrix();

        for (Size i=0; i<3; ++i)
            for (Size j=0; j<3; ++j)
                EXPECT_NEAR(FLOAT_INIT, sigmamat(i,j), ZERO_EPS);

    };


}
