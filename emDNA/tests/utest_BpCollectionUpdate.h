// BpCollectionUpdate unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_BpCollectionUpdate_h
#define emDNA_utest_BpCollectionUpdate_h


#include <gtest/gtest.h>


namespace {


	// testing class for BpCollectionUpdate
	class BpCollectionUpdateTest : public ::testing::Test {

	protected:

		BpCollectionUpdateTest() {};
		virtual ~BpCollectionUpdateTest() {};

		// set up before test
		// called after the constructor
		virtual void SetUp() {

            // all results have been computed in Mathematica

            // bp set 1
            Triad bp1;
            bp1.set_origin(Vector3(0,0,0));
            bp1.set_axis(I, Vector3(1,0,0));
            bp1.set_axis(J, Vector3(0,1,0));
            bp1.set_axis(K, Vector3(0,0,1));
            Triad bp2;
            bp2.set_origin(Vector3(0.3619202045165672,-0.5587498313793305,
                                   2.4866558872098694));
            bp2.set_axis(I, Vector3(0.8464643144004624,0.5225525835217476,
                                    0.10216145017227315));
            bp2.set_axis(J, Vector3(-0.5283277263480836,0.8481287954579576,
                                    0.03933649561028737));
            bp2.set_axis(K, Vector3(-0.06609068026900065,-0.08727166647761882,
                                    0.9939897777199723));
            Triad bp3;
            bp3.set_origin(Vector3(-1.4695350696613938,-0.2055433224796379,
                                   4.344999022927246));
            bp3.set_axis(I, Vector3(0.4444746120822953,0.8955582886300123,
                                    0.02043699783179022));
            bp3.set_axis(J, Vector3(-0.8946029677107045,0.4449442010677555,
                                    -0.04135442055423358));
            bp3.set_axis(K, Vector3(-0.046128617771324584,
                                    0.00009799112231430318,0.9989355039341864));
            Triad bp4;
            bp4.set_origin(Vector3(-3.0773733748317245,1.0053124051842495,
                                   8.543033635461985));
            bp4.set_axis(I, Vector3(-0.22003738644555518,0.9738957029770319,
                                    -0.05577372400227519));
            bp4.set_axis(J, Vector3(-0.9695590919648419,-0.22463745070891236,
                                    -0.09743296632717154));
            bp4.set_axis(K, Vector3(-0.10741841441075177,0.03263702593487504,
                                    0.993678070998654));
            Triad bp5;
            bp5.set_origin(Vector3(-2.635650729958132,1.673368563131532,
                                   11.450855188057718));
            bp5.set_axis(I, Vector3(-0.6703097623187417,0.7412992572887493,
                                    -0.034062203148647906));
            bp5.set_axis(J, Vector3(-0.7418487369839653,-0.6705426524761077,
                                    0.005744792911962601));
            bp5.set_axis(K, Vector3(-0.01858154932955887,0.029119793156103774,
                                    0.999403203752649));
            bp_set_1.push_back(bp1);
            bp_set_1.push_back(bp2);
            bp_set_1.push_back(bp3);
            bp_set_1.push_back(bp4);
            bp_set_1.push_back(bp5);

            // bp set 2
            Triad optbp1;
            optbp1.set_origin(Vector3(0,0,0));
            optbp1.set_axis(I, Vector3(1,0,0));
            optbp1.set_axis(J, Vector3(0,1,0));
            optbp1.set_axis(K, Vector3(0,0,1));
            Triad optbp2;
            optbp2.set_origin(Vector3(7.800752408471633E-8,7.40148352746332E-8,
                                      3.400000008378934));
            optbp2.set_axis(I, Vector3(0.8260982938612407,0.5635260498677477,
                                       1.2132158631439494E-9));
            optbp2.set_axis(J, Vector3(-0.5635260498677477,0.8260982938612407,
                                       -7.832681771719644E-9));
            optbp2.set_axis(K, Vector3(-5.416155773316892E-9,
                                       5.786886304981238E-9,
                                       1.));
            Triad optbp3;
            optbp3.set_origin(Vector3(5.518369526723006E-7,5.620950931407218E-7,
                                      6.800000005512624));
            optbp3.set_axis(I, Vector3(0.3648767828324339,0.9310558164524041,
                                       4.327327930662956E-9));
            optbp3.set_axis(J, Vector3(-0.931055816452404,0.36487678283243397,
                                       -2.15075528611554E-8));
            optbp3.set_axis(K, Vector3(-2.1603673682637515E-8,
                                       3.818622855036204E-9,
                                       0.9999999999999998));
            Triad optbp4;
            optbp4.set_origin(Vector3(-1.3746496237900925E-8,
                                      1.0948831472667991E-6,
                                      10.199999966656877));
            optbp4.set_axis(I, Vector3(-0.22325011720913276,0.9747611939168016,
                                       7.788819940083086E-9));
            optbp4.set_axis(J, Vector3(-0.9747611939168009,-0.22325011720913204,
                                       -4.369439612800079E-8));
            optbp4.set_axis(K, Vector3(-4.0852746772659345E-8,
                                       -1.7347018480956817E-8,
                                       0.999999999999999));
            Triad optbp5;
            optbp5.set_origin(Vector3(-5.7101167352200164E-8,
                                      1.3933178402703452E-7,
                                      13.599999959503105));
            optbp5.set_axis(I, Vector3(-0.7337298656262844,0.6794413030483464,
                                       -1.494596946099064E-8));
            optbp5.set_axis(J, Vector3(-0.6794413030483446,-0.7337298656262835,
                                       -5.897379778203019E-8));
            optbp5.set_axis(K, Vector3(-5.103553817499947E-8,
                                       -3.311592775618442E-8,
                                       0.9999999999999981));
            bp_set_2.push_back(optbp1);
            bp_set_2.push_back(optbp2);
            bp_set_2.push_back(optbp3);
            bp_set_2.push_back(optbp4);
            bp_set_2.push_back(optbp5);

            // updated dofs
            Real udofs1[6], udofs2[6], udofs3[6], udofs4[6];
            udofs1[0]=-7.126663638835746E-9;
            udofs1[1]=-3.468928054386503E-9;
            udofs1[2]=0.5986479345600425;
            udofs1[3]=7.800752408471633E-8;
            udofs1[4]=7.40148352746332E-8;
            udofs1[5]=3.400000008378934;
            udofs2[0]=-1.1433068361879199E-8;
            udofs2[1]=-1.1627327504872407E-8;
            udofs2[2]=0.5986479339247115;
            udofs2[3]=4.7382942858758426E-7;
            udofs2[4]=4.880802578660885E-7;
            udofs2[5]=3.39999999713369;
            udofs3[0]=-1.762758801065415E-8;
            udofs3[1]=-2.2533959459348057E-8;
            udofs3[2]=0.598647933046616;
            udofs3[3]=-5.655834489102015E-7;
            udofs3[4]=5.327880541260773E-7;
            udofs3[5]=3.3999999611442533;
            udofs4[0]=-1.671047859290072E-8;
            udofs4[1]=-8.550300911087194E-9;
            udofs4[2]=0.5986479338582795;
            udofs4[3]=-4.335467111429924E-8;
            udofs4[4]=-9.555513632397646E-7;
            udofs4[5]=3.3999999928462272;
            updated_dofs.push_back(BpStepDofs(VectorN(Array<Real>(6,
                                                                  udofs1))));
            updated_dofs.push_back(BpStepDofs(VectorN(Array<Real>(6,
                                                                  udofs2))));
            updated_dofs.push_back(BpStepDofs(VectorN(Array<Real>(6,
                                                                  udofs3))));
            updated_dofs.push_back(BpStepDofs(VectorN(Array<Real>(6,
                                                                  udofs4))));

            // update params
            Real uprms1[6], uprms2[6], uprms3[6], uprms4[6];
            uprms1[0]=-4.0832774851463386E-7;
            uprms1[1]=-1.9875493695087473E-7;
            uprms1[2]=34.30000006451433;
            uprms1[3]=1.0226125208273072E-7;
            uprms1[4]=3.5606101187121116E-8;
            uprms1[5]=3.400000008378934;
            uprms2[0]=-6.550665640202278E-7;
            uprms2[1]=-6.661967930455671E-7;
            uprms2[2]=34.30000002811254;
            uprms2[3]=6.928546912803244E-7;
            uprms2[4]=-1.1249445069126448E-7;
            uprms2[5]=3.399999997133686;
            uprms3[0]=-1.0099863960058937E-6;
            uprms3[1]=-1.291100772739542E-6;
            uprms3[2]=34.29999997780138;
            uprms3[3]=5.202120187867675E-7;
            uprms3[4]=4.993327344460325E-7;
            uprms3[5]=3.3999999611442657;
            uprms4[0]=-9.57439897016922E-7;
            uprms4[1]=-4.898961557721587E-7;
            uprms4[2]=34.30000002430627;
            uprms4[3]=-8.093712084701406E-7;
            uprms4[4]=3.378529621687832E-7;
            uprms4[5]=3.399999992846249;
            updated_prms.push_back(BpStepParams(VectorN(Array<Real>(6,
                                                                    uprms1))));
            updated_prms.push_back(BpStepParams(VectorN(Array<Real>(6,
                                                                    uprms2))));
            updated_prms.push_back(BpStepParams(VectorN(Array<Real>(6,
                                                                    uprms3))));
            updated_prms.push_back(BpStepParams(VectorN(Array<Real>(6,
                                                                    uprms4))));

            
        };

		// clean up before test
		// called after the destructor
		virtual void TearDown() {};

		// declare objects that will be used by all tests
        std::vector<BasePair> bp_set_1;
        std::vector<BasePair> bp_set_2;
        std::vector<BpStepDofs> updated_dofs;
        std::vector<BpStepParams> updated_prms;
        
	};
    
    
}


#endif  // emDNA_utest_BpCollectionUpdate_h
