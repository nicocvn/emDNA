// BpCollection unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_BpCollection_h
#define emDNA_utest_BpCollection_h


#include <gtest/gtest.h>


namespace {


	// testing class for BpCollection
	class BpCollectionTest : public ::testing::Test {

	protected:

		BpCollectionTest() {};
		virtual ~BpCollectionTest() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {

            // all results have been computed in Mathematica

            // bp step parameters
            Real p1[6], p2[6], p3[6];
            p1[0]=-3.381968273358222;
            p1[1]=4.5688010062250655;
            p1[2]=37.350487852227474;
            p1[3]=-1.4148348990327229;
            p1[4]=-1.286198360073577;
            p1[5]=1.8663760624888797;
            p2[0]=-0.8869414726289051;
            p2[1]=-2.8848660034386064;
            p2[2]=30.806545866466696;
            p2[3]=0.8101391712679584;
            p2[4]=-0.909199483811777;
            p2[5]=4.150147935390688;
            p3[0]=-0.038955502318522406;
            p3[1]=-6.926156900732028;
            p3[2]=41.773392099156396;
            p3[3]=-1.4854452499214372;
            p3[4]=-0.14555886721685773;
            p3[5]=2.697919924811686;
            test_bp_step_params.assign(3, BpStepParams());
            test_bp_step_params[0] =
            (BpStepParams(VectorN(Array<Real>(6, p1))));
            test_bp_step_params[1] =
            (BpStepParams(VectorN(Array<Real>(6, p2))));
            test_bp_step_params[2] =
            (BpStepParams(VectorN(Array<Real>(6, p3))));

            // bp step dofs - numerical values
            Real d1[6], d2[6], d3[6];
            d1[0]=-0.05902648156808837;
            d1[1]=0.07974062042705701;
            d1[2]=0.6518889902475147;
            d1[3]=-0.8743248588064237;
            d1[4]=-1.593613702062678;
            d1[5]=1.9584113534902727;
            d2[0]=-0.015480048969861077;
            d2[1]=-0.050350410238852124;
            d2[2]=0.5376756565364929;
            d2[3]=1.359282448730724;
            d2[4]=0.36758843122974527;
            d2[5]=4.089489825591035;
            d3[0]=-0.000679901777226536;
            d3[1]=-0.12088424242749873;
            d3[2]=0.7290832318568641;
            d3[3]=0.15518770802387438;
            d3[4]=-1.5136972437589409;
            d3[5]=2.681630907591572;
            test_bp_step_dofs.assign(3, BpStepDofs());
            test_bp_step_dofs[0] = (BpStepDofs(VectorN(Array<Real>(6, d1))));
            test_bp_step_dofs[1] = (BpStepDofs(VectorN(Array<Real>(6, d2))));
            test_bp_step_dofs[2] = (BpStepDofs(VectorN(Array<Real>(6, d3))));

            // base pairs
            Triad bp1;
            bp1.set_origin(Vector3(0,0,0));
            bp1.set_axis(I, Vector3(1,0,0));
            bp1.set_axis(J, Vector3(0,1,0));
            bp1.set_axis(K, Vector3(0,0,1));
            Triad bp2;
            bp2.set_origin(Vector3(-0.8743248588064237,-1.593613702062678,
                                   1.9584113534902727));
            bp2.set_axis(I, Vector3(0.7922666851501041,0.6028460057246604,
                                    -0.09428781990845193));
            bp2.set_axis(J, Vector3(-0.6075489546142423,0.7937027272034514,
                                    -0.030335599169050077));
            bp2.set_axis(K, Vector3(0.05654880501308,0.08131835101394147,
                                    0.9950826892474646));
            Triad bp3;
            bp3.set_origin(Vector3(0.48495758992430027,-1.2260252708329327,
                                   6.0479011790813075));
            bp3.set_axis(I, Vector3(0.3708921214646523,0.927206637822158,
                                    -0.052219584582600564));
            bp3.set_axis(J, Vector3(-0.9285729346911856,0.37110177127940897,
                                    -0.005981664675049521));
            bp3.set_axis(K, Vector3(0.013832541142145736,0.05070824520543927,
                                    0.9986177119767784));
            Triad bp4;
            bp4.set_origin(Vector3(0.6401452979481747,-2.7397225145918735,
                                   8.72953208667288));
            bp4.set_axis(I, Vector3(-0.3405848845717334,0.9376277539285744,
                                    0.06968594882840992));
            bp4.set_axis(J, Vector3(-0.940147195706596,-0.3405060756361166,
                                    -0.013373962012621025));
            bp4.set_axis(K, Vector3(0.011188690999521167,-0.0700700186795182,
                                    0.9974793259391241));
            test_base_pairs.assign(4, BasePair());
            test_base_pairs[0] = bp1;
            test_base_pairs[1] = bp2;
            test_base_pairs[2] = bp3;
            test_base_pairs[3] = bp4;

        };

		// clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests
        std::vector<BpStepDofs> test_bp_step_dofs;
        std::vector<BpStepParams> test_bp_step_params;
        std::vector<BasePair> test_base_pairs;
        
	};
    
    
}


#endif  // emDNA_utest_BpCollection_h
