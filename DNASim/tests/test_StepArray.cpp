// test_StepArray class
// Nicolas Clauvelin


#include "test_StepArray.h"
#include <DNASim.h>
using namespace DNASim;


namespace {
	
	
	// CopyOperator
	TEST_F(StepArrayTest, CopyOperator) {
		
		StepParameters p;
		p[TILT]=Real(2.35);
		p[ROLL]=Real(0.78);
		p[TWIST]=Real(32.52);
		p[SHIFT]=Real(-0.87); 
		p[SLIDE]=Real(-0.61);
		p[RISE]=Real(3.07);
		
		StepArray steps;
		steps.append_step_parameters(p);
		steps.append_step_parameters(p);
		steps.append_step_parameters(p);
		
		// copy
		StepArray steps_bis = steps;
		
		// checking
		EXPECT_EQ(3, (Integer)steps_bis.n_step());
		EXPECT_EQ(4, (Integer)steps_bis.n_bp());
		EXPECT_EQ(steps.n_step(), steps_bis.n_step());
		EXPECT_EQ(steps.n_bp(), steps_bis.n_bp());
		EXPECT_EQ(steps.parameters(2)[SHIFT], steps_bis.parameters(2)[SHIFT]);
		EXPECT_EQ(steps.parameters(2)[TILT], steps_bis.parameters(2)[TILT]);
		
	};
	
	
	// MethodBasePair
	TEST_F(StepArrayTest, MethodBasePair) {
		
		StepParameters p;
		p[TILT]=Real(2.35);
		p[ROLL]=Real(0.78);
		p[TWIST]=Real(32.52);
		p[SHIFT]=Real(-0.87); 
		p[SLIDE]=Real(-0.61);
		p[RISE]=Real(3.07);
		
		StepArray steps;
		steps.append_step_parameters(p);
		steps.append_step_parameters(p);
		steps.append_step_parameters(p);
        steps.append_step_parameters(p);
		
		// before last bp
		Triad lastbp = steps.base_pair(3, Triad());
        Triad lastbp_bis = steps.base_pairs(0, 3, Triad()).back();

		// pre-computed result
		Triad res;
		res.set_origin(Vector3(Real(0.1260121690678413),
                               Real(-3.1728193734702774),
                               Real(9.05380378885885)));
		res.set_axis(I, Vector3(Real(-0.13025529453219345),
                                Real(0.9897394784626935),
                                Real(0.058730937492276925)));
		res.set_axis(J, Vector3(Real(-0.9857156740761646),
                                Real(-0.1356490144274099),
                                Real(0.09981961112652465)));
		res.set_axis(K, Vector3(Real(0.1067622036439408),
                                Real(-0.044889972791950145),
                                Real(0.9932707194998888)));
		
		// checking
        EXPECT_GE(ZERO_TOL, (res.origin()-lastbp.origin()).norm());
        EXPECT_GE(ZERO_TOL, (res.axis(I)-lastbp.axis(I)).norm());
        EXPECT_GE(ZERO_TOL, (res.axis(J)-lastbp.axis(J)).norm());
        EXPECT_GE(ZERO_TOL, (res.axis(K)-lastbp.axis(K)).norm());
        EXPECT_GE(ZERO_TOL, (res.origin()-lastbp_bis.origin()).norm());
        EXPECT_GE(ZERO_TOL, (res.axis(I)-lastbp_bis.axis(I)).norm());
        EXPECT_GE(ZERO_TOL, (res.axis(J)-lastbp_bis.axis(J)).norm());
        EXPECT_GE(ZERO_TOL, (res.axis(K)-lastbp_bis.axis(K)).norm());

        // check full reconstruction size
        EXPECT_TRUE(steps.all_base_pairs(Triad()).size()==5);

	};
	
	
}

