// AlglibMinimizer unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_AlglibMinimizer_h
#define emDNA_utest_AlglibMinimizer_h


#include <gtest/gtest.h>


namespace {


	// testing class for AlglibMinimizer
	class AlglibMinimizerTest : public ::testing::Test {

	protected:

		AlglibMinimizerTest() {};
		virtual ~AlglibMinimizerTest() {};

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
            bp2.set_origin(Vector3(0.3423522534064986,0.7194586291377235,
                                   2.8641802332928794));
            bp2.set_axis(I, Vector3(0.8936729534272733,0.28539738425273065,
                                    -0.3462614407846678));
            bp2.set_axis(J, Vector3(-0.3229608283649402,0.9448252003250356,
                                    -0.05478726286818889));
            bp2.set_axis(K, Vector3(0.3115203936412614,0.16079077676426476,
                                    0.9365369028784386));
            Triad bp3;
            bp3.set_origin(Vector3(2.1305689366365064,2.571458123608726,
                                   6.1709717247426905));
            bp3.set_axis(I, Vector3(0.5461161449243888,0.5370863045052877,
                                    -0.6428805937075543));
            bp3.set_axis(J, Vector3(-0.7662707054000991,0.6303818830527974,
                                    -0.12428953119407826));
            bp3.set_axis(K, Vector3(0.33850607424174634,0.5604970856584918,
                                    0.7558152252169696));
            Triad bp4;
            bp4.set_origin(Vector3(2.4836653314871846,4.54561390457075,
                                   9.90045375624199));
            bp4.set_axis(I, Vector3(0.21952766866315312,0.8318237861671759,
                                    -0.5097811211273184));
            bp4.set_axis(J, Vector3(-0.9158497380865949,0.35577305212192395,
                                    0.1861311167713944));
            bp4.set_axis(K, Vector3(0.3361946756539083,0.426021976135452,
                                    0.8399276254009078));
            Triad bp5;
            bp5.set_origin(Vector3(4.230749844621096,6.203510448988075,
                                   12.225244146779392));
            bp5.set_axis(I, Vector3(-0.4244465856470388,0.8523431105957789,
                                    -0.3055426611006768));
            bp5.set_axis(J, Vector3(-0.7809695746276278,-0.17386228771688866,
                                    0.599882012078869));
            bp5.set_axis(K, Vector3(0.4581829541116882,0.49323739388834065,
                                    0.7394492909129907));
            test_base_pairs.assign(5, BasePair());
            test_base_pairs[0] = bp1;
            test_base_pairs[1] = bp2;
            test_base_pairs[2] = bp3;
            test_base_pairs[3] = bp4;
            test_base_pairs[4] = bp5;

            // imposed last bp gradients
            Real impg1[6], impg2[6], impg3[6];
            impg1[0]=28.61367209582602;
            impg1[1]=40.58189391862218;
            impg1[2]=-88.6653367456207;
            impg1[3]=-11.749520471094918;
            impg1[4]=7.046190743806502;
            impg1[5]=-2.2267360867736894;
            impg2[0]=-17.87346467867345;
            impg2[1]=46.259975867369576;
            impg2[2]=-50.77397578580875;
            impg2[3]=5.276968910675083;
            impg2[4]=10.114511408897613;
            impg2[5]=14.556261881391038;
            impg3[0]=17.54686272254494;
            impg3[1]=-31.57216065484373;
            impg3[2]=-47.676379356829436;
            impg3[3]=-23.78311518156443;
            impg3[4]=4.048224024433337;
            impg3[5]=27.565005229856254;
            mm_imp_gradients.assign(3, VectorN(6, FLOAT_INIT));
            mm_imp_gradients[0] = VectorN(Array<Real>(6, impg1));
            mm_imp_gradients[1] = VectorN(Array<Real>(6, impg2));
            mm_imp_gradients[2] = VectorN(Array<Real>(6, impg3));

        };

        // clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests
        std::vector<BasePair> test_base_pairs;
        std::vector<VectorN> mm_imp_gradients;

	};
    
    
}


#endif  // emDNA_utest_AlglibMinimizer_h
