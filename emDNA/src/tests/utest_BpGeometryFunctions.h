// BpGeometricFunctions unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_BpGeometryFunctions_h
#define emDNA_utest_BpGeometryFunctions_h


#include <gtest/gtest.h>
#include <BpGeometryFunctions.h>
using namespace BpGeometryFunctions;


namespace {


	// testing class for BpGeometryFunctions
	class BpGeometryFunctionsTest : public ::testing::Test {

	protected:

		BpGeometryFunctionsTest() {};
		virtual ~BpGeometryFunctionsTest() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {

            // all results have been computed in Mathematica

            // ThetaAngles and ZYZEulerAngles structures
            test_theta_angles._Theta1Degree = Real(-11.0);
            test_theta_angles._Theta2Degree = Real(7.56);
            test_theta_angles._Theta3Degree = Real(42.15);
            test_theta_angles_degen._Theta1Degree = FLOAT_INIT;
            test_theta_angles_degen._Theta2Degree = FLOAT_INIT;
            test_theta_angles_degen._Theta3Degree = Real(38.26);
            test_euler_angles._ZetaRadian = Real(1.3364914653228088);
            test_euler_angles._KappaRadian = Real(0.23295641214122734);
            test_euler_angles._EtaRadian = Real(-0.6008351856071991);

            // Xi matrix
            Real xielems[9];
            xielems[0] = 0.9319243612494583;
            xielems[1] = -0.3613155946356177;
            xielems[2] = 0.026798454058496048;
            xielems[3] = 0.35462020106400444;
            xielems[4] = 0.925879777518865;
            xielems[5] = 0.11227359486237368;
            xielems[6] = -0.06567562721277274;
            xielems[7] = -0.09555977504503994;
            xielems[8] = 0.9864940726840916;
            mm_ximat << Array<Real>(9, xielems);

            // Xi matrix limit
            Real xielems_limit[9];
            xielems_limit[0] = 0.9447774517425483;
            xielems_limit[1] = -0.3277126281955836;
            xielems_limit[2] = FLOAT_INIT;
            xielems_limit[3] = 0.3277126281955836;
            xielems_limit[4] = 0.9447774517425483;
            xielems_limit[5] = FLOAT_INIT;
            xielems_limit[6] = FLOAT_INIT;
            xielems_limit[7] = FLOAT_INIT;
            xielems_limit[8] = Real(1);
            mm_ximat_limit << Array<Real>(9, xielems_limit);

            // Lambda vectors
            mm_Lambda1 = Vector3(-0.4996375380048692,
                                 0.0005273917918569901,
                                 -0.016474722505744735);
            mm_Lambda2 = Vector3(0.0005273917918569901,
                                 -0.49923263099068427,
                                 -0.023971157085078312);
            mm_Lambda3 = Vector3(0.03291218401824212,
                                 0.04788809843924117,
                                 -0.49661203989736585);

            // dLambda matrices
            Real dl1_elems[9], dl2_elems[9], dl3_elems[9];

            dl1_elems[0] = 0.000001739748859;
            dl1_elems[1] = 0.005492863195;
            dl1_elems[2] = Real(0.0);
            dl1_elems[3] = -0.002744498058;
            dl1_elems[4] = 0.003995260757;
            dl1_elems[5] = Real(0.0);
            dl1_elems[6] = -0.0001318181288;
            dl1_elems[7] = -0.1247681442;
            dl1_elems[8] = Real(0.0);

            dl2_elems[0] = -0.002744498058;
            dl2_elems[1] = 0.003995260757;
            dl2_elems[2] = Real(0.0);
            dl2_elems[3] = -0.007990317786;
            dl2_elems[4] = -0.000002531380615;
            dl2_elems[5] = Real(0.0);
            dl2_elems[6] = 0.1246669403;
            dl2_elems[7] = 0.0001318181288;
            dl2_elems[8] = Real(0.0);

            dl3_elems[0] = 0.0005270340195;
            dl3_elems[1] = 0.2490728684;
            dl3_elems[2] = Real(0.0);
            dl3_elems[3] = -0.2486682361;
            dl3_elems[4] = -0.0005270340195;
            dl3_elems[5] = Real(0.0);
            dl3_elems[6] = -0.02394404922;
            dl3_elems[7] = 0.01645609201;
            dl3_elems[8] = Real(0.0);

            mm_dL1 << Array<Real>(9, dl1_elems);
            mm_dL2 << Array<Real>(9, dl2_elems);
            mm_dL3 << Array<Real>(9, dl3_elems);

            // dXi matrices
            Real dXi1_elems[9], dXi2_elems[9], dXi3_elems[9];

            dXi1_elems[0] = -0.00790913649770305;
            dXi1_elems[1] = -0.0025183026106408114;
            dXi1_elems[2] = 0.1798942646336981;
            dXi1_elems[3] = 0.02036914416269704;
            dXi1_elems[4] = 0.06730048638794744;
            dXi1_elems[5] = -0.4551356832163941;
            dXi1_elems[6] = -0.0021033724031344053;
            dXi1_elems[7] = 0.4946824285359707;
            dXi1_elems[8] = 0.09512722501086969;

            dXi2_elems[0] = -0.029438611082711796;
            dXi2_elems[1] = -0.029753706267272126;
            dXi2_elems[2] = 0.4611638632495397;
            dXi2_elems[3] = -0.04548362423149324;
            dXi2_elems[4] = -0.011379443932765363;
            dXi2_elems[5] = 0.17321704091717124;
            dXi2_elems[6] = -0.49629730053933124;
            dXi2_elems[7] = 0.0021033724031344053;
            dXi2_elems[8] = -0.06537834737110682;
            
            dXi3_elems[0] = -0.17731010053200227;
            dXi3_elems[1] = -0.4629398887594325;
            dXi3_elems[2] = -0.05613679743118684;
            dXi3_elems[3] = 0.46596218062472916;
            dXi3_elems[4] = -0.18065779731780884;
            dXi3_elems[5] = 0.013399227029248024;
            dXi3_elems[6] = 0;
            dXi3_elems[7] = 0;
            dXi3_elems[8] = 0;

            mm_dXi1 << Array<Real>(9, dXi1_elems);
            mm_dXi2 << Array<Real>(9, dXi2_elems);
            mm_dXi3 << Array<Real>(9, dXi3_elems);

        };

		// clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests
        ThetaAngles test_theta_angles;
        ThetaAngles test_theta_angles_degen;
        ZYZEulerAngles test_euler_angles;
        Matrix3 mm_ximat, mm_ximat_limit;
        Vector3 mm_Lambda1, mm_Lambda2, mm_Lambda3;
        Matrix3 mm_dL1, mm_dL2, mm_dL3;
        Matrix3 mm_dXi1, mm_dXi2, mm_dXi3;

	};


}


#endif  // emDNA_utest_BpGeometryFunctions_h
