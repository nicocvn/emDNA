// BpGeometry functions
// Nicolas Clauvelin


// step geometry functions implemented in a sub-namespace

// not designed to be used directly, see StepParameters class
//
// see supporting notes for implementation details


#ifndef DNASim_BpGeometry_h
#define DNASim_BpGeometry_h


#include <DNASim_Includes.h>


namespace DNASim {
	
	
	class Vector3;
	class Matrix4;
	class StepParameters;
	
	
namespace BpGeometry {	
	
	
	// step frame computation from bp triads function
	Triad step_frame_from_triads(const Triad& bp1, const Triad& bp2);
	
	// step parameters computation from bp triads function
	StepParameters bp_step_parameters_from_triads(const Triad& bp1,
												  const Triad& bp2);
	
	// bp triad reconstruction function
	Triad reconstruct_bp_from_parameters(const Triad& bp,
										 const StepParameters& p);
	
	// step generator matrix computation from parameters function
	Matrix4 step_matrix_from_parameters(const StepParameters& p);
	
	// local step displacement from parameters function
	Vector3 local_step_displacement(const StepParameters& p);
	
	// Euler angles computation from step parameters function
	Vector3 euler_zyz_angles_from_parameters(const StepParameters& p);
	
	
}


}


#endif	// DNASim_BpGeometry_h
