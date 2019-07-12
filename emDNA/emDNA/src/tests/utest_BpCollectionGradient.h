// BpCollectionGradient unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_BpCollectionGradient_h
#define emDNA_utest_BpCollectionGradient_h


#include <gtest/gtest.h>


namespace {


	// testing class for BpCollectionGradient
	class BpCollectionGradientTest : public ::testing::Test {

	protected:

		BpCollectionGradientTest() {};
		virtual ~BpCollectionGradientTest() {};

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

            // full gradients
            Real g1[6], g2[6], g3[6], g4[6];
            g1[0]=0.038622294111583244;
            g1[1]=-36.53057423613445;
            g1[2]=-56.01455367567452;
            g1[3]=-2.144639524919092;
            g1[4]=1.1084245360335927;
            g1[5]=16.09764317792072;
            g2[0]=-64.64275597174444;
            g2[1]=7.900988376193881;
            g2[2]=-25.659189574166632;
            g2[3]=0.4614881387656995;
            g2[4]=3.5438141469435105;
            g2[5]=10.043072174389412;
            g3[0]=-61.37759281570407;
            g3[1]=-118.60247959287621;
            g3[2]=37.38809930091578;
            g3[3]=-12.105929210363772;
            g3[4]=3.339270798886647;
            g3[5]=14.232489931540115;
            g4[0]=37.55741106757519;
            g4[1]=117.9181595740985;
            g4[2]=-37.38076058774386;
            g4[3]=-7.116410756720558;
            g4[4]=1.6558826384311707;
            g4[5]=12.46221018853573;
            mm_gradients.assign(4, VectorN(6, FLOAT_INIT));
            mm_gradients[0] = VectorN(Array<Real>(6, g1));
            mm_gradients[1] = VectorN(Array<Real>(6, g2));
            mm_gradients[2] = VectorN(Array<Real>(6, g3));
            mm_gradients[3] = VectorN(Array<Real>(6, g4));

            // imposed last bp gradients
            Real impg1[6], impg2[6], impg3[6];
            impg1[0]=84.42698690749847;
            impg1[1]=-40.06793076710525;
            impg1[2]=42.01757803830219;
            impg1[3]=4.971771231801465;
            impg1[4]=-0.5474581023975709;
            impg1[5]=3.6354329893849986;
            impg2[0]=14.780883902735752;
            impg2[1]=-45.318025807512136;
            impg2[2]=60.271250467946004;
            impg2[3]=7.577898895486263;
            impg2[4]=1.8879315085123411;
            impg2[5]=-2.419138014146327;
            impg3[0]=-27.092987949755546;
            impg3[1]=-229.24736648895788;
            impg3[2]=92.352678129379;
            impg3[3]=-4.989518453643198;
            impg3[4]=1.6833881604554826;
            impg3[5]=1.7702797430043962;
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
        std::vector<VectorN> mm_gradients, mm_imp_gradients;
        
	};
    
    
}


#endif  // emDNA_utest_BpCollectionGradient_h
