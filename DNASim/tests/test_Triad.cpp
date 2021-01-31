// test_Triad class
// Nicolas Clauvelin


#include "test_Triad.h"
#include <DNASim.h>
using namespace DNASim;


namespace {
	
	
	// MethodOrthogonalize
	TEST_F(TriadTest, MethodOrthogonalize) {
		
		Triad f;
		f.set_axis(I, Vector3(Real(0.98),Real(0.1),Real(-0.06)));
		f.set_axis(J, Vector3(Real(0.04),Real(0.96),Real(0.12)));
		f.set_axis(K, Vector3(Real(-0.2),Real(-0.16),Real(1.16)));
		
		f.orthogonalize();

        // axis
        Vector3 axis_I = f.axis(I);
        Vector3 axis_J = f.axis(J);
        Vector3 axis_K = f.axis(K);

        // results
        Vector3 resI(Real(0.9929939380407246),
                     Real(0.1013259120449719),
                     Real(-0.06079554722698314));
        Vector3 resJ(Real(-0.09252926340439611),
                     Real(0.9867505536779633),
                     Real(0.13327295385813573));
		Vector3 resK(Real(0.07349404348798344),
                     Real(-0.12671386808273005),
                     Real(0.9892129301658459));
		
		EXPECT_GE(ZERO_EPS, (resI-axis_I).norm());
        EXPECT_GE(ZERO_EPS, (resJ-axis_J).norm());
        EXPECT_GE(ZERO_EPS, (resK-axis_K).norm());
		
	};
	
	
	// MethodRotate
	TEST_F(TriadTest, MethodRotate) {

        // rotation transformation
        Real angle = Real(0.25);
        Vector3 axis(Real(2.10),Real(-1.03),Real(3.078));
        axis.normalize();
        AffineTransformation t = AffineTransformation::rotation(angle, axis);

        // tiriad
		Triad f2;
		f2.set_axis(I, Vector3(Real(0.843184),Real(0.537622),
                               Real(-0.00158442)));
		f2.set_axis(J, Vector3(Real(-0.537064),Real(0.842436),
                               Real(0.043173)));
		f2.set_axis(K, Vector3(Real(0.0245455),Real(-0.0355518),
                               Real(0.999066)));
		f2.transform(t);
		
		// pre-compute result
		Vector3 x(Real(0.7164684086555715),
						 Real(0.6846171364277268),
						 Real(0.13405811812341092));
		Vector3 y(Real(-0.6972955261028181),
						 Real(0.7086431672875318),
						 Real(0.10772127391880294));
		Vector3 z(Real(-0.021251576606688638),
						 Real(-0.17065697750662562),
						 Real(0.9851009343866867));
        
		// axis
        Vector3 axis_I = f2.axis(I);
        Vector3 axis_J = f2.axis(J);
        Vector3 axis_K = f2.axis(K);

        EXPECT_GE(ZERO_EPS, (axis_I-x).norm());
        EXPECT_GE(ZERO_EPS, (axis_J-y).norm());
        EXPECT_GE(ZERO_EPS, (axis_K-z).norm());
		
	};
	
	
	// EulerZYZAngles
	TEST_F(TriadTest, EulerZYZAngles) {
		
		Triad f1;
		Triad f2;
		f2.set_origin(Vector3(Real(-0.626757), Real(-0.883724), Real(3.0627)));
		f2.set_axis(I, Vector3(Real(0.843184), Real(0.537622),
                               Real(-0.00158442)));
		f2.set_axis(J, Vector3(Real(-0.537064), Real(0.842436),
                               Real(0.043173)));
		f2.set_axis(K, Vector3(Real(0.0245455), Real(-0.0355518),
                               Real(0.999066)));
		
		Vector3 eangles = Triad::euler_zyz_angles(f1, f2);
		
		// pre-computed result
		Vector3 res(Real(-0.9665321805647735), Real(0.04322373026220664),
                    Real(1.5341134635722116));

        EXPECT_GE(ZERO_EPS, (res-eangles).norm());

        // full transformation checking
        f2.orthogonalize();
        f1.transform(Triad::euler_zyz_transformation(f1, f2));

        // difference matrix
        Matrix4 m = f1.matrix_representation()-f2.matrix_representation();

		// check all entries
		for (Size i=0; i<3; ++i)
			for (Size j=0; j<3; ++j)
				EXPECT_GE(ZERO_EPS, std::sqrt(m(i,j)*m(i,j)));

        // zero angles case
        f1.orthogonalize();
        eangles = Triad::euler_zyz_angles(f1, f1);
        for (Size i=0; i<3; ++i)
            EXPECT_GE(ZERO_TOL_SQ,
                      eangles[(CartesianCoordinates)i]
                      *eangles[(CartesianCoordinates)i]);

	};


    // EulerZYZAnglesTorture
	TEST_F(TriadTest, EulerZYZAnglesTorture) {

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
        const Vector3 zero_angles = Triad::euler_zyz_angles(f1, f1);
        for (Size i=0; i<3; ++i)
            EXPECT_NEAR(zero_angles[(CartesianCoordinates)i],
                        FLOAT_INIT,
                        ZERO_EPS);

        // modified zero case
        Triad f2(f1);
        f2.transform(AffineTransformation::
                     rotation(ZERO_EPS, Vector3(1.0, FLOAT_INIT, FLOAT_INIT)));
        const Vector3 zero_angles_mod = Triad::euler_zyz_angles(f1, f2);
        for (Size i=0; i<3; ++i)
            EXPECT_NEAR(zero_angles_mod[(CartesianCoordinates)i],
                        FLOAT_INIT,
                        ZERO_EPS);

        // flipped case
        Triad f_trivial;
        Triad f_flipped;
        f_flipped.transform(AffineTransformation::
                            rotation(F_PI/Real(2.0),
                                     Vector3(1.0, FLOAT_INIT, FLOAT_INIT)));
        const Vector3 flipped_angles = Triad::euler_zyz_angles(f_trivial,
                                                               f_flipped);
        EXPECT_NEAR(flipped_angles[X]+flipped_angles[Z],
                    FLOAT_INIT,
                    ZERO_EPS);
        EXPECT_NEAR(flipped_angles[Y]-flipped_angles[Z],
                    FLOAT_INIT,
                    ZERO_EPS);
        EXPECT_NEAR(flipped_angles[Y]+flipped_angles[X],
                    FLOAT_INIT,
                    ZERO_EPS);

	};


	// TransformationMatrix
	TEST_F(TriadTest, TransformationMatrix) {

        Triad f1;
        Triad f2;

        // f1
        f1.set_origin(Vector3(Real(-1.3201), Real(-4.26086), Real(19.832)));
        f1.set_axis(I,
                    Vector3(Real(-0.732544), Real(-0.650647), Real(-0.200095)));
        f1.set_axis(J,
                    Vector3(Real(0.664275), Real(-0.747487), Real(-0.0013067)));
        f1.set_axis(K,
                    Vector3(Real(-0.148718), Real(-0.133875), Real(0.979776)));

        // f2
        f2.set_origin(Vector3(Real(-22.6723), Real(-36.3599), Real(74.7274)));
        f2.set_axis(I,
                    Vector3(Real(-0.909892), Real(0.358993), Real(-0.207895)));
        f2.set_axis(J,
                    Vector3(Real(-0.21324), Real(-0.834602), Real(-0.507905)));
        f2.set_axis(K,
                    Vector3(Real(-0.355844), Real(-0.417807), Real(0.83595)));

        // transformation matrix
        Matrix4 tmat = Triad::transformation_matrix(f1, f2);
        Matrix4 tmatinv = Triad::transformation_matrix(f2, f1);

        // first checking
        Matrix4 f2mat = f2.matrix_representation();
        Matrix4 f2m = f1.matrix_representation()*tmat;
        for (Size i=0; i<4; ++i)
            for (Size j=0; j<4; ++j)
                EXPECT_GE(ZERO_TOL, std::fabs(f2mat(i,j)-f2m(i,j)));

        // second checking
        Matrix4 f1mat = f1.matrix_representation();
        Matrix4 f1m = f2mat*tmatinv;
        for (Size i=0; i<4; ++i)
            for (Size j=0; j<4; ++j)
                EXPECT_GE(ZERO_TOL, std::fabs(f1mat(i,j)-f1m(i,j)));
        
	};


    // PointTransformation
    TEST_F(TriadTest, PointTransformation) {

        // frame 1
        Triad f1;

        // frame 2 (rotation of angle pi/2)
        Triad f2;
        f2.set_origin(Vector3(0, 1, 0));
        f2.set_axis(I, Vector3(1/std::sqrt(2), 1/std::sqrt(2), 0));
        f2.set_axis(J, Vector3(-1/std::sqrt(2), 1/std::sqrt(2), 0));

        // point
        Vector3 pt(1,0,0);
        Vector3 pt_in_f2 = f2.express_point_in(pt);
        Vector3 pt_back_in_f1 = f2.express_point_out(pt_in_f2);

        // checking
        EXPECT_GE(ZERO_EPS, (pt_back_in_f1-pt).norm());
        EXPECT_GE(ZERO_EPS, (pt_in_f2-Vector3(0,-2/std::sqrt(2),0)).norm());

    };


    // VectorTransformation
    TEST_F(TriadTest, VectorTransformation) {

        // frame 1
        Triad f1;

        // frame 2 (rotation of angle pi/2)
        Triad f2;
        f2.set_origin(Vector3(0, 1, 0));
        f2.set_axis(I, Vector3(1/std::sqrt(2), 1/std::sqrt(2), 0));
        f2.set_axis(J, Vector3(-1/std::sqrt(2), 1/std::sqrt(2), 0));

        // point
        Vector3 pt(1,0,0);
        Vector3 pt_in_f2 = f2.express_vector_in(pt);
        Vector3 pt_back_in_f1 = f2.express_vector_out(pt_in_f2);

        // checking
        EXPECT_GE(ZERO_EPS, (pt_back_in_f1-pt).norm());
        EXPECT_GE(ZERO_EPS, (pt_in_f2-Vector3(1/std::sqrt(2),
                                              -1/std::sqrt(2),0)).norm());

    };
	
	
	
}

