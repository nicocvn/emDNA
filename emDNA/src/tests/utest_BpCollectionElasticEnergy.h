// BpCollectionElasticEnergy unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_BpCollectionElasticEnergy_h
#define emDNA_utest_BpCollectionElasticEnergy_h


#include <gtest/gtest.h>


namespace {


	// testing class for BpCollectionElasticEnergy
	class BpCollectionElasticEnergyTest : public ::testing::Test {

	protected:

		BpCollectionElasticEnergyTest() {};
		virtual ~BpCollectionElasticEnergyTest() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {

            // all results have been computed in Mathematica

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

            // precomputed energy
            exactE = 61.716460378123635;

        };

		// clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests
        std::vector<BasePair> test_base_pairs;
        Real exactE;
        
	};
    
    
}


#endif  // emDNA_utest_BpCollectionElasticEnergy_h
