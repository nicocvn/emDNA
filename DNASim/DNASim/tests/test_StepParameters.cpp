// test_StepParameters class
// Nicolas Clauvelin


#include <Vector3.h>
#include <AffineTransformation.h>
#include <Triad.h>
#include <BpGeometry.h>
#include <StepParameters.h>
#include <test_StepParameters.h>
using namespace DNASim;


namespace {
	
	
	// MethodStepLocalDisplacement
	TEST_F(StepParametersTest, MethodStepLocalDisplacement) {

        StepParameters p;
		p[TILT]=Real(2.35);
		p[ROLL]=Real(0.78);
		p[TWIST]=Real(32.52);
		p[SHIFT]=Real(-0.87); 
		p[SLIDE]=Real(-0.61);
		p[RISE]=Real(3.07);
		
		Vector3 rvec = p.step_local_displacement();
		
		// pre-computed results
		Vector3 res(Real(-0.6267572231993119),
						   Real(-0.8837238877554395), 
						   Real(3.0626961118233553));
		
        EXPECT_GE(ZERO_EPS, (res-rvec).norm());
		
	};
	
	
	// MethodStepMatrix
	TEST_F(StepParametersTest, MethodStepMatrix) {
		
		StepParameters p;
		p[TILT]=Real(2.35);
		p[ROLL]=Real(0.78);
		p[TWIST]=Real(32.52);
		p[SHIFT]=Real(-0.87); 
		p[SLIDE]=Real(-0.61);
		p[RISE]=Real(3.07);
		
		Matrix4 M = p.step_matrix();
		
		// pre-computed result
		Real vals[] = {
			Real(0.843184387199365), Real(-0.537063875605327),
			Real(0.02454552307483492), Real(-0.6267572231993118),
			Real(0.5376221524421768), Real(0.842436043444874),
			Real(-0.03555184817118858), Real(-0.8837238877554396),
			Real(-0.0015844199796980927), Real(0.04317298026234038),
			Real(0.999066355848597), Real(3.0626961118233553),
			Real(0.0), Real(0.0), Real(0.0), Real(1.0)
		};
		Matrix4 res;
        res << Array<Real>(16, vals);
		
		// check all entries
		for (Size i=0; i<4; ++i)
			for (Size j=0; j<4; ++j)
				EXPECT_GE(ZERO_EPS, fabs(res(i,j)-M(i,j)));
		
	};
	
	
	// MethodStepEulerZYZAngles
	TEST_F(StepParametersTest, MethodStepEulerZYZAngles) {
		
		StepParameters p;
		p[TILT]=Real(2.35);
		p[ROLL]=Real(0.78);
		p[TWIST]=Real(32.52);
		p[SHIFT]=Real(-0.87); 
		p[SLIDE]=Real(-0.61);
		p[RISE]=Real(3.07);
		
		Vector3 angles = p.step_euler_zyz_angles();
		
		// pre-computed result
		Vector3 res(Real(-0.9665323745377767),
						   Real(0.04321549420023254),
						   Real(1.534113447286333));

        EXPECT_GE(ZERO_EPS, (res-angles).norm());
		
	};


    // degenerated case
	TEST_F(StepParametersTest, DegeneratedEulerZYZAngles) {

		StepParameters p;
		p[TILT]=Real(0.0); p[ROLL]=Real(0.0); p[TWIST]=Real(45.0);
		p[SHIFT]=Real(0.0); p[SLIDE]=Real(0.0); p[RISE]=Real(3.4);

        // rebuild bp and check axis
        Triad bp = p.reconstruct_triad(Triad());
        EXPECT_GE(ZERO_TOL, std::fabs((bp.axis(I))[X]-Real(1)/std::sqrt(2)));
        EXPECT_GE(ZERO_TOL, std::fabs((bp.axis(I))[Y]-Real(1)/std::sqrt(2)));

        StepParameters pbis(Triad(), bp);
        EXPECT_GE(ZERO_TOL, std::fabs(pbis[TWIST]-Real(45.0)));
        EXPECT_GE(ZERO_TOL, std::fabs(pbis[RISE]-Real(3.4)));

	};


    // StepParametersTorture
	TEST_F(StepParametersTest, StepParametersTorture) {

		// zero case
        Triad f1;
		f1.set_origin(Vector3(Real(-0.626757), Real(-0.883724), Real(3.0627)));
		f1.set_axis(I, Vector3(Real(0.843184), Real(0.537622),
                               Real(-0.00158442)));
		f1.set_axis(J, Vector3(Real(-0.537064), Real(0.842436),
                               Real(0.043173)));
		f1.set_axis(K, Vector3(Real(0.0245455), Real(-0.0355518),
                               Real(0.999066)));
        f1.orthogonalize();
        const StepParameters p0 = StepParameters(f1, f1);
        for (Size i=0; i<3; ++i)
            EXPECT_NEAR(p0[(BpStepParameters)i],
                        FLOAT_INIT,
                        ZERO_TOL);

        // modified zero case
        Triad f2(f1);
        f2.transform(AffineTransformation::
                     rotation(ZERO_EPS, Vector3(1.0, FLOAT_INIT, FLOAT_INIT)));
        const StepParameters p0_mod = StepParameters(f1, f2);
        for (Size i=0; i<3; ++i)
            EXPECT_NEAR(p0_mod[(BpStepParameters)i],
                        FLOAT_INIT,
                        ZERO_TOL);

        // flipped case
        Triad f_trivial;
        Triad f_flipped;
        f_flipped.transform(AffineTransformation::
                            rotation(F_PI/Real(2.0),
                                     Vector3(FLOAT_INIT, 1.0, FLOAT_INIT)));
        const StepParameters p_flipped = StepParameters(f_trivial,
                                                        f_flipped);

	};

	
}



