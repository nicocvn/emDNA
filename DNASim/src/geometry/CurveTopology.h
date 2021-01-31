// CurveTopology class
// Nicolas Clauvelin


// in this class a curve is always assumed to be closed, that is, the last and
// first point always form the last segment of the curve

// implementation of methods and formula described in:
// Characterization of the Geometry and Topology of DNA Pictured As a Discrete
// Collection of Atoms
// N. Clauvelin, W. K. Olson, I. Tobias
// JCTC, 2012, DOI: 10.1021/ct200657e


#ifndef DNASim_CurveTopology_h
#define DNASim_CurveTopology_h


#include "geometry/Segment.h"


namespace DNASim {
	
	
	typedef Vector3 Vertex;
	typedef std::vector<Vertex> CurvePoints;
	
	
	namespace {
		class SegmentQuadruplet {
		public:
			SegmentQuadruplet() : s_imm(), s_im(), s_i(), s_ip() {};
			Segment s_imm, s_im, s_i, s_ip;
		};
	}
	
	
	class CurveTopology {
		
		
	public:
		
		// curve segments computation methods
        // the curve is always assumed to be closed
		static std::vector<Segment> curve_segments(const CurvePoints& pts);
		static Segment				curve_segment(const CurvePoints& pts,
												  Size i);
		
		// curve writhing number
		// pts is the list of points describing the curve
		// curve is closed by joining the last and first point
		static Real curve_writhing_number(const CurvePoints& pts);
		
		// curve writhing number gradient
		// pts is the list of points describing the curve
		// the calculation is done assumign the curve is closed;
		// if not discard values for the first two and last two points
		static
        std::vector<Vector3> curve_writhing_gradient(const CurvePoints& pts);
		
		// two curves linking number
		// pts1 is the list of points describing the first curve
		// pts2 is the list of points describing the second curve
		// curves are closed by joining the last and first point
		static Real two_curves_linking_number(const CurvePoints& pts1,
											  const CurvePoints& pts2);


        // dihedral angle computation function
		static Real dihedral_angle(Real s1, const Segment& segment1,
								   Real s2, const Segment& segment2);

		
	private:
		

		// vertex writhe gradient 
		static Vector3 writhe_gradient(const SegmentQuadruplet& sq);

        // deleted constructors and operators
        CurveTopology() = delete;
        ~CurveTopology() = delete;
        CurveTopology(const CurveTopology& curve_topo) = delete;
        CurveTopology(CurveTopology&& curve_topo) = delete;
        CurveTopology& operator=(const CurveTopology& curve_topo) = delete;
        CurveTopology& operator=(CurveTopology&& curve_topo) = delete;
		
		
	};
	
	
}


#endif	// DNASim_CurveTopology_h
