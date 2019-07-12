// CurveTopology class
// Nicolas Clauvelin


#include <CurveTopology.h>
#include <Segment.h>


namespace DNASim {
	
	
	// curve segments computation methods
    // the curve is always assumed to be closed
	std::vector<Segment>
	CurveTopology::curve_segments(const CurvePoints& pts) {
		
		Size n_pts = pts.size();
		std::vector<Segment> segments;
		segments.reserve(n_pts);
		
		// curve segments
		for (Size i=0; i<n_pts-1; ++i)
			segments.push_back(Segment(pts[i],pts[i+1]));
		
		// closing segment
		segments.push_back(Segment(pts[n_pts-1],pts[0]));
		
		return segments;
		
	};
	Segment CurveTopology::curve_segment(const CurvePoints& pts, Size i) {
		
		if (i != pts.size()-1)
			return Segment(pts[i],pts[i+1]);
		else
			return Segment(pts[pts.size()-1],pts[0]);
		
	};
	
	
	// curve writhing number
	Real CurveTopology::curve_writhing_number(const CurvePoints& pts) {
		
		Real Wr = FLOAT_INIT;
		
		// number of vertices
		Size n = pts.size();
		
		// curve segments
		const std::vector<Segment> s = CurveTopology::curve_segments(pts);
		
		// loop over segments
		for (Size i=0; i<n; ++i)
			for (Size j=i+2; j<n; ++j) {
				Wr += CurveTopology::dihedral_angle(Real(1.0), s[i],
													Real(0.0), s[j]);
				Wr += CurveTopology::dihedral_angle(Real(0.0), s[i],
													Real(1.0), s[j]);
				Wr += -CurveTopology::dihedral_angle(Real(1.0), s[i],
													 Real(1.0), s[j]);
				Wr += -CurveTopology::dihedral_angle(Real(0.0), s[i],
													 Real(0.0), s[j]);
			};
		
		return Wr/(Real(2.0)*F_PI);
		
	};
	
	
	// curve writhing number gradient
	// pts is the list of points describing the curve
	// the calculation is done assuming the curve is closed;
	// if not discard values for the first two and last two points
	std::vector<Vector3>
	CurveTopology::curve_writhing_gradient(const CurvePoints& pts) {
		
		std::vector<Vector3> gradients;
		gradients.reserve(pts.size());
		
		// curve segments
		Size n_pts = pts.size();
		std::vector<Segment> segments = curve_segments(pts);
		
		// create quadruplets
		std::vector<SegmentQuadruplet> squads(n_pts, SegmentQuadruplet());
		std::rotate(segments.begin(), segments.begin()+(Integer)n_pts-2,
                    segments.end());
		for (Size i=0; i<n_pts; ++i)
			squads[i].s_imm = segments[i];
		std::rotate(segments.begin(), segments.begin()+1, segments.end());
		for (Size i=0; i<n_pts; ++i)
			squads[i].s_im = segments[i];
		std::rotate(segments.begin(), segments.begin()+1, segments.end());
		for (Size i=0; i<n_pts; ++i)
			squads[i].s_i = segments[i];
		std::rotate(segments.begin(), segments.begin()+1, segments.end());
		for (Size i=0; i<n_pts; ++i)
			squads[i].s_ip = segments[i];
		
		// computation
		for (auto it = squads.begin(), end = squads.end(); it != end; ++it)
			gradients.push_back(writhe_gradient(*it));
		
		return gradients;
		
	};
	
	
	// two curves linking number
	Real CurveTopology::two_curves_linking_number(const CurvePoints& pts1,
												  const CurvePoints& pts2) {
		
		Real Lk = FLOAT_INIT;
		
		// loop over pairs of segments
        for (const Segment& s_i : CurveTopology::curve_segments(pts1))
            for (const Segment& s_j : CurveTopology::curve_segments(pts2)) {
                Lk += CurveTopology::dihedral_angle(Real(1.0), s_i,
													Real(0.0), s_j);
				Lk += CurveTopology::dihedral_angle(Real(0.0), s_i,
													Real(1.0), s_j);
				Lk += -CurveTopology::dihedral_angle(Real(1.0), s_i,
													 Real(1.0), s_j);
				Lk += -CurveTopology::dihedral_angle(Real(0.0), s_i,
													 Real(0.0), s_j);
            };

		return Lk/(Real(4.0)*F_PI);
		
	};
	
		
	// dihedral angle computation function
	Real CurveTopology::dihedral_angle(Real s1, const Segment& segment1,
									   Real s2, const Segment& segment2) {
	
		// relative lengths
		const Real L1 = segment1.length()*s1;
		const Real L2 = segment2.length()*s2;
		
		// segment tangents
		const Vector3 t1 = segment1.tangent();
		const Vector3 t2 = segment2.tangent();
		
		// joining vector for the two points
		const Vector3 rjoin = (segment1.first()+L1*t1)-(segment2.first()+L2*t2);
		
		// dihedral norm
		const Real norm1 = (rjoin.cross(t1)).norm();
		const Real norm2 = (rjoin.cross(t2)).norm();
		const Real norm = norm1*norm2;
		
		// dihedral angle
		if (norm1*norm1 < ZERO_EPS || norm2*norm2 < ZERO_EPS)
			return FLOAT_INIT;
		else {
		
			// cosine of the dihedral angle
			Real cosine = rjoin.cross(t2).dot(t1.cross(rjoin))/(norm*norm);
			
			// sine of the dihedral angle
			Real sine = Real(-1.0)*(rjoin.norm())*(rjoin.dot(t1.cross(t2)))
                /(norm*norm);

			return std::atan2(sine, cosine);
			
		};
		
	};
	
	
	// vertex writhe gradient 
	Vector3 CurveTopology::writhe_gradient(const SegmentQuadruplet& s) {
		return Real(-1.0/(2.0*F_PI))*
			(
			 Real(1.0/s.s_im.length())*Segment::binormal_vector(s.s_imm,s.s_im)
			 *std::tan(Segment::turning_angle(s.s_imm,s.s_im)*Real(0.5))
			 +
			 (Real(1.0/s.s_im.length())*Segment::binormal_vector(s.s_im,s.s_i)
			  -
			  Real(1.0/s.s_i.length())*Segment::binormal_vector(s.s_im,s.s_i))
			 *std::tan(Segment::turning_angle(s.s_im,s.s_i)*Real(0.5))
			 +
			 Real(-1.0/s.s_i.length())*Segment::binormal_vector(s.s_i,s.s_ip)
			 *std::tan(Segment::turning_angle(s.s_i,s.s_ip)*Real(0.5))
			 );
	};
	
	
}

