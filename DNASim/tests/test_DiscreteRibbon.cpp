// test_DiscreteRibon class
// Nicolas Clauvelin


#include <Segment.h>
#include <CurveTopology.h>
#include <Triad.h>
#include <DiscreteRibbon.h>
#include <test_DiscreteRibbon.h>
using namespace DNASim;


namespace {
	
	
	// RibbonTotalTwistTest
	TEST_F(DiscreteRibbonTest, RibbonTopologyTest) {
		
		std::vector<Vertex> basecurve(20, Vertex());
        basecurve[0] = {0.0, 0.0, 0.0};
        basecurve[1] = {0.5058452589888371, -0.16435908789150877,
            3.3581403580234683};
        basecurve[2] = {1.9738653775743604, -0.6413477388861121,
            6.387562540263919};
        basecurve[3] = {4.2603603182518155, -1.3842749804491077,
            8.79172559629818};
        basecurve[4] = {7.141512025280681, -2.3204179178978457,
            10.33529329541264};
        basecurve[5] = {10.335293295412642, -3.3581403580234674,
            10.867170476549425};
        basecurve[6] = {13.529074565544601, -4.395862798149089,
            10.335293295412642};
        basecurve[7] = {16.410226272573468, -5.3320057355978285,
            8.791725596298184};
        basecurve[8] = {18.696721213250925, -6.074932977160826,
            6.387562540263923};
        basecurve[9] = {20.16474133183645, -6.55192162815543,
            3.3581403580234737};
        basecurve[10] = {20.670586590825287, -6.71628071604694, 0};
        basecurve[11] = {20.16474133183645, -6.551921628155433,
            -3.3581403580234626};
        basecurve[12] = {18.696721213250925, -6.074932977160831,
            -6.387562540263913};
        basecurve[13] = {16.41022627257347, -5.3320057355978365,
            -8.791725596298175};
        basecurve[14] = {13.529074565544605, -4.395862798149099,
            -10.335293295412635};
        basecurve[15] = {10.335293295412644, -3.358140358023477,
            -10.867170476549422};
        basecurve[16] = {7.141512025280683, -2.320417917897855,
            -10.335293295412638};
        basecurve[17] = {4.260360318251816, -1.3842749804491161,
            -8.791725596298182};
        basecurve[18] = {1.9738653775743606, -0.6413477388861192,
            -6.3875625402639225};
        basecurve[19] = {0.5058452589888367, -0.16435908789151432,
            -3.358140358023474};

		std::vector<Vertex> edgecurve(20, Vertex());
        edgecurve[0] = {0., 2., 0.};
        edgecurve[1] = {-0.5944089685709038, 1.4292031590059628,
            3.8581403580234683};
        edgecurve[2] = {0.4350236089867334, -0.14134773888611196,
            7.563133044848865};
        edgecurve[3] = {2.992580557061661, -2.2084163430527655,
            10.100742590673129};
        edgecurve[4] = {6.372091140986868, -4.070417917897846,
            10.923078547705114};
        edgecurve[5] = {9.747508043120169, -5.167157352398415,
            10.24913648779953};
        edgecurve[6] = {12.690275043394344, -5.359388289711511,
            8.796451526825015};
        edgecurve[7] = {15.292192283823573, -4.968734471595149,
            7.173691607548289};
        edgecurve[8] = {17.815043334812216, -4.552390491223457,
            5.436506023968769};
        edgecurve[9] = {20.193509589753976, -4.5612690017791415,
            3.1671573523984202};
        edgecurve[10] = {21.846157095410234, -5.098246727297044,
            0.};
        edgecurve[11] = {21.99153808740155, -5.909415897553114,
            -3.858140358023462};
        edgecurve[12] = {20.23556298183855, -6.574932977160831,
            -7.563133044848858};
        edgecurve[13] = {16.951463505758266, -6.743932350493969,
            -10.100742590673121};
        edgecurve[14] = {13.122924945253471, -6.263896786898995,
            -10.923078547705106};
        edgecurve[15] = {9.74750804312017, -5.167157352398425,
            -10.249136487799527};
        edgecurve[16] = {7.253769019425579, -3.5929604038352236,
            -8.796451526825011};
        edgecurve[17] = {5.378394307001711, -1.7475462444517964,
            -7.173691607548287};
        edgecurve[18] = {3.5820857840184317, 0.07217775267630222,
            -5.436506023968769};
        edgecurve[19] = {1.6526475056562582, 1.4630222744820935,
            -3.167157352398422};

		// total twist checking
		DiscreteRibbon ribbon(basecurve, edgecurve);
		Real TW = ribbon.ribbon_total_twist();
		Real mm_TW = Real(2.0);
		EXPECT_LE(std::fabs(TW-mm_TW), ZERO_TOL);
		
		// twist + writhe
		// linking number decomposition precomputed to 2
		Real WR = ribbon.ribbon_writhing_number();
		EXPECT_LE(std::fabs(TW+WR-2.0), ZERO_TOL);
		
		// direct linking number precomputed to 2
		Real LK = ribbon.ribbon_linking_number();
		EXPECT_LE(std::fabs(LK-2.0), ZERO_TOL);
		
	};
	
	
}
