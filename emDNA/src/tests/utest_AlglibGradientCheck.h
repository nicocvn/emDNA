// AlglibGradientCheck unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_AlglibGradientCheck_h
#define emDNA_utest_AlglibGradientCheck_h


#include <gtest/gtest.h>


namespace {


	// testing class for AlglibGradientCheck
	class AlglibGradientCheckTest : public ::testing::Test {

	protected:

		AlglibGradientCheckTest() {};
		virtual ~AlglibGradientCheckTest() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {

            // base pairs
            Triad bp1;
            bp1.set_origin(Vector3(0,0,0));
            bp1.set_axis(I, Vector3(1,0,0));
            bp1.set_axis(J, Vector3(0,1,0));
            bp1.set_axis(K, Vector3(0,0,1));
            Triad bp2;
            bp2.set_origin(Vector3(-0.5722699393814894,0.2060196169533935,
                                   4.067820738857911));
            bp2.set_axis(I, Vector3(0.9492245014825983,0.30169276487671337,
                                    0.08918700245078874));
            bp2.set_axis(J, Vector3(-0.2992119522251905,0.9533349288896988,
                                    -0.040307828079212786));
            bp2.set_axis(K, Vector3(-0.09718566473870101,0.011575360897923687,
                                    0.9951989537722493));
            Triad bp3;
            bp3.set_origin(Vector3(-0.5714539118243419,1.337147571045705,
                                   7.827080225574836));
            bp3.set_axis(I, Vector3(0.7349719905459917,0.6768517980120091,
                                    -0.041083045660870926));
            bp3.set_axis(J, Vector3(-0.6342840908309417,0.6647963617152167,
                                    -0.39462575887665624));
            bp3.set_axis(K, Vector3(-0.2397912951539912,0.3160972017878702,
                                    0.917923032606901));
            Triad bp4;
            bp4.set_origin(Vector3(-3.2724258567824265,1.935265685247943,
                                   11.124672153695935));
            bp4.set_axis(I, Vector3(-0.0024684960425339214,0.9878179787896915,
                                    -0.15559417504243853));
            bp4.set_axis(J, Vector3(-0.7884672248348426,-0.09762480507482696,
                                    -0.6072798636504786));
            bp4.set_axis(K, Vector3(-0.6150718184801917,0.12118183945604655,
                                    0.7791030868232516));
            Triad bp5;
            bp5.set_origin(Vector3(-5.336692942448856,2.3790622639577,
                                   14.636686885396454));
            bp5.set_axis(I, Vector3(-0.16443339959895112,0.7837724027359432,
                                    -0.5988843609628307));
            bp5.set_axis(J, Vector3(-0.8627205358887562,-0.4086207460844021,
                                    -0.29789656397018727));
            bp5.set_axis(K, Vector3(-0.47819968010460756,0.46768569208278443,
                                    0.743367445728489));
            test_base_pairs.assign(5, BasePair());
            test_base_pairs[0] = bp1;
            test_base_pairs[1] = bp2;
            test_base_pairs[2] = bp3;
            test_base_pairs[3] = bp4;
            test_base_pairs[4] = bp5;

        };

		// clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests
        std::vector<BasePair> test_base_pairs;

	};


}


#endif  // emDNA_utest_AlglibGradientCheck_h
