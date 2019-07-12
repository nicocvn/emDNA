// BpGeometry functions
// Nicolas Clauvelin


#include <Vector3.h>
#include <Matrix4.h>
#include <AffineTransformation.h>
#include <Triad.h>
#include <StepParameters.h>
#include <BpGeometry.h>


namespace DNASim {
namespace BpGeometry {
	
	
	// step frame computation from bp triads function
	Triad step_frame_from_triads(const Triad& bp1, const Triad& bp2) {
		
		// euler angles
		const Vector3 eangles(Triad::euler_zyz_angles(bp1, bp2));
		
		// angles for the mid step frame
		const Real gamma = (eangles[Z]-eangles[X])/Real(2.0);
		const Vector3 midangles(eangles[X], eangles[Y]/Real(2.0), gamma);

        // midstep frame
        Triad midstep_frame(bp1);
        midstep_frame.transform(Triad::euler_zyz_transformation(bp1,
                                                                midangles));
		
		return midstep_frame;

	};
	
	
	// step parameters computation from bp triads functions
	StepParameters bp_step_parameters_from_triads(const Triad& bp1,
												  const Triad& bp2) {
		
		StepParameters parameters;
		
		// euler angles
		const Vector3 eangles(Triad::euler_zyz_angles(bp1, bp2));
		
		// gamma angle
		const Real gamma = (eangles[Z]-eangles[X])/Real(2.0);
		
		// rotation bp step parameters
		parameters[TILT] = eangles[Y]*std::sin(gamma)*RAD_2_DEG;
		parameters[ROLL] = eangles[Y]*std::cos(gamma)*RAD_2_DEG;
		parameters[TWIST] = (eangles[X]+eangles[Z])*RAD_2_DEG;
		
		// joining vector
		const Vector3 rvec(bp2.origin()-bp1.origin());
		
		// mid step frame
		const Triad midframe(step_frame_from_triads(bp1, bp2));
		
		// transport in mid step frame
		const Vector3 rho(midframe.express_vector_in(rvec));
		
		// translation bp step parameters
		parameters[SHIFT] = rho[X];
		parameters[SLIDE] = rho[Y];
		parameters[RISE] = rho[Z];
		
		return parameters;
		
	};
	
	
	// bp triad reconstruction function
	Triad reconstruct_bp_from_parameters(const Triad& bp,
										 const StepParameters& p) {
		
		// euler angles
		const Vector3 eangles(euler_zyz_angles_from_parameters(p));
		
		// rotation
        Triad frame(bp);
		frame.transform(Triad::euler_zyz_transformation(bp, eangles));
		
		// translation vector
		frame.set_origin(bp.origin()
                         +bp.express_vector_out(local_step_displacement(p)));
		
		return frame;
		
	};
	
	
	// step generator matrix computation from parameters function
	Matrix4 step_matrix_from_parameters(const StepParameters& p) {
		
		// euler angles computation
		const Vector3 eangles(euler_zyz_angles_from_parameters(p));
		
		// step rotation
        Triad rotframe;
		rotframe.transform(Triad::euler_zyz_transformation(Triad(), eangles));
		Matrix4 Grot(rotframe.matrix_representation());
		
		// mid-step frame
		const Vector3 midangles(eangles[X], eangles[Y]/Real(2.0),
                                (eangles[Z]-eangles[X])/Real(2.0));
        Triad midframe;
		midframe.transform(Triad::euler_zyz_transformation(Triad(), midangles));
		
		// step translation
		const Vector3 Gtrans(p[SHIFT]*midframe.axis(I)
                             + p[SLIDE]*midframe.axis(J)
                             + p[RISE]*midframe.axis(K));
		
		// generator matrix construction
		const Real vals[] = {
			Grot(0,0), Grot(0,1), Grot(0,2), Gtrans[X],
			Grot(1,0), Grot(1,1), Grot(1,2), Gtrans[Y],
			Grot(2,0), Grot(2,1), Grot(2,2), Gtrans[Z],
			FLOAT_INIT, FLOAT_INIT, FLOAT_INIT, 1.0
		};
		Matrix4 G;
        G << Array<Real>(16, vals);
		
		return G;

	};
	
	
	// local step displacement from parameters function
	Vector3 local_step_displacement(const StepParameters& p) {
		
		// euler angles
		const Vector3 eangles(euler_zyz_angles_from_parameters(p));
		
		// angles for the mid step frame
		const Real gamma = (eangles[Z]-eangles[X])/Real(2.0);
		const Vector3 midangles(eangles[X], eangles[Y]/Real(2.0), gamma);
		
		// mid step frame
        Triad midframe;
		midframe.transform(Triad::euler_zyz_transformation(Triad(), midangles));
		
		return p[SHIFT]*midframe.axis(I)
            + p[SLIDE]*midframe.axis(J) + p[RISE]*midframe.axis(K);

	};
	
	
	// Euler angles computation from step parameters function
	Vector3 euler_zyz_angles_from_parameters(const StepParameters& p) {
		
		// rotational bp step parameters
		const Real vtilt = p[TILT]*DEG_2_RAD;
		const Real vroll = p[ROLL]*DEG_2_RAD;
		const Real vtwist = p[TWIST]*DEG_2_RAD;
        const Real gamma = std::atan2(vtilt, vroll);
		
		// euler angles
        return Vector3(vtwist/Real(2.0) - gamma,
                       std::hypot(vtilt, vroll),
                       vtwist/Real(2.0) + gamma);

	};
	
	
}
}

