// test_CurveTopology class
// Nicolas Clauvelin


#include <CurveTopology.h>
#include <test_CurveTopology.h>
using namespace DNASim;


namespace {
	
	
	// LinkingNumberTest
	TEST_F(CurveTopologyTest, LinkingNumberTest) {
		
		std::vector<Vertex> pts1;
		std::vector<Vertex> pts2;
		
		pts1.push_back(Vector3(0.0,0.0,0.0));
		pts1.push_back(Vector3(1.0,0.0,0.0));
		pts1.push_back(Vector3(1.0,1.0,0.0));
		pts1.push_back(Vector3(0.0,1.0,0.0));
		
		pts2.push_back(Vector3(0.0,0.0,1.0));
		pts2.push_back(Vector3(1.0,0.0,1.0));
		pts2.push_back(Vector3(1.0,1.0,1.0));
		pts2.push_back(Vector3(0.0,1.0,1.0));
		
		Real Lk0 = CurveTopology::two_curves_linking_number(pts1, pts2);
		
		// linking number is zero
		EXPECT_TRUE(Lk0*Lk0<=ZERO_TOL_SQ);

		pts2.push_back(Vector3(Real(-0.5),Real(0.2),Real(-0.5)));
		pts2.push_back(Vector3(Real(0.5),Real(0.2),Real(-0.5)));
		pts2.push_back(Vector3(Real(0.5),Real(0.2),Real(0.5)));
		pts2.push_back(Vector3(Real(-0.5),Real(0.2),Real(0.5)));
		
		Real Lk1 = CurveTopology::two_curves_linking_number(pts1, pts2);
		
		// linking number is one
		EXPECT_TRUE((Lk1-Real(1.0))*(Lk1-Real(1.0))<=ZERO_TOL_SQ);
		
	};
	
	
	// WrithingNumberTest
	TEST_F(CurveTopologyTest, WrithingNumberTest) {
		
		// arbitrary curve vertices
        std::vector<Vertex> curve(10, Vertex());
        curve[0] = {4.690960256616805, -9.530717643860346, 1.0887294887904275};
        curve[1] = {4.922627441332729, -0.6643562857377816, 8.832933209905555};
        curve[2] = {-6.013939770113872, 5.876999938404282, 6.698559385686529};
        curve[3] = {-2.446909444815091, -4.03576511187698, -1.95648576184972};
        curve[4] = {-2.1086163901059436, -8.238200094876568,
            -9.501808606743523};
        curve[5] = {0.7745176479225648, 1.0180267803201133, 4.335655829889916};
        curve[6] = {-0.3767259537912935, 8.675102728925154, 2.3444101031478226};
        curve[7] = {-6.23396606574531, -2.9436018169086644,
            -1.6448827929040668};
        curve[8] = {-0.6038392252251796, -3.963801534645132,
            -8.557313894380506};
        curve[9] = {6.013589944831942, -2.9923659939655813, 0.9362585344378331};

		// Mathematica result
		Real WrTest = Real(-0.6290795961316287);

		Real Wr = CurveTopology::curve_writhing_number(curve);
		EXPECT_TRUE(std::fabs(WrTest-Wr)<ZERO_EPS);
		
	};
	
	
}


