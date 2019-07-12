// DiscreteRibbon class
// Nicolas Clauvelin


#include <Segment.h>
#include <CurveTopology.h>
#include <Triad.h>
#include <AffineTransformation.h>
#include <DiscreteRibbon.h>


namespace DNASim {
	
	
	// class default constructor
	DiscreteRibbon::DiscreteRibbon(const CurvePoints& base_curve,
								   const CurvePoints& edge_reference_curve) :
	m_ribbon_curves(base_curve, edge_reference_curve) {};
	
	
	// class destructor
	DiscreteRibbon::~DiscreteRibbon() {};
	
	
	// ribbon curves accessors
	const CurvePoints& DiscreteRibbon::ribbon_base_curve() const {
		return m_ribbon_curves.first;
	};
	const CurvePoints& DiscreteRibbon::ribbon_edge_reference_curve() const {
		return m_ribbon_curves.second;	
	};
	
	
	// ribbon topology methods
	Real DiscreteRibbon::ribbon_writhing_number() const {
		return CurveTopology::curve_writhing_number(ribbon_base_curve());
	};
	Real DiscreteRibbon::ribbon_linking_number() const {
		return CurveTopology::
				two_curves_linking_number(ribbon_base_curve(),
										  ribbon_edge_reference_curve());
	};
	Real DiscreteRibbon::ribbon_total_twist() const {
		
		// ribbon material frames and segments
		std::vector<Triad> mframes = ribbon_material_frames();
		std::vector<Segment> segments =
			CurveTopology::curve_segments(ribbon_base_curve());
		
		// cyclization
		mframes.insert(mframes.begin(), mframes.back());
		segments.insert(segments.begin(), segments.back());
		
		// twist computation
		Real Tw = FLOAT_INIT;
		for (Size i=0; i<ribbon_base_curve().size(); ++i)
			Tw += vertex_twist(i)*vertex_length(i);
		
		return Tw/(Real(2.0)*F_PI);
		
	};
	
	
	// writhe gradient method
	std::vector<Vector3> DiscreteRibbon::ribbon_writhe_gradient() const {
		return CurveTopology::curve_writhing_gradient(m_ribbon_curves.first);
	};
	
	
	// ribbon vertex property
	Real DiscreteRibbon::vertex_curvature(Size i) const {
		const Segment sim =
        CurveTopology::curve_segment(ribbon_base_curve(), i-1);
		const Segment si =
        CurveTopology::curve_segment(ribbon_base_curve(), i);
		return Segment::turning_angle(sim, si)/vertex_length(i);
	};
	Real DiscreteRibbon::vertex_twist(Size i) const {
		
		// incoming and outgoing frames
		Triad mat1;
		if (i != 0)
			mat1 = ribbon_material_frame(i-1);
		else
			mat1 = ribbon_material_frame(ribbon_base_curve().size()-1);
		Triad mat2 = ribbon_material_frame(i);
		
		// incoming and outgoing segments
		Segment sim;
		if (i != 0)
			sim = CurveTopology::curve_segment(ribbon_base_curve(), i-1);
		else
			sim = CurveTopology::curve_segment(ribbon_base_curve(),
											   ribbon_base_curve().size()-1);
		Segment si = CurveTopology::curve_segment(ribbon_base_curve(), i);
		
		// parallel transport rotation
        mat1.transform(Segment::parallel_transport_rotation(sim, si));

		// twist computation
        return std::atan2(mat1.axis(I).cross(mat2.axis(I)).dot(mat2.axis(K)),
                          mat1.axis(I).dot(mat2.axis(I)))/vertex_length(i);

	};
	Real DiscreteRibbon::vertex_length(Size i) const {
		Segment sim;
		if (i != 0)
			sim = CurveTopology::curve_segment(ribbon_base_curve(), i-1);
		else
			sim = CurveTopology::curve_segment(ribbon_base_curve(),
											   ribbon_base_curve().size()-1);
		Segment si = CurveTopology::curve_segment(ribbon_base_curve(), i);
		return Real(0.5)*(sim.length()+si.length());
	};
	
	
	// ribbon all material frames accessor
	std::vector<Triad> DiscreteRibbon::ribbon_material_frames() const {
		
		// material frames building
		std::vector<Triad> material_frames;
		material_frames.reserve(m_ribbon_curves.first.size());
		for (Size i=0; i<m_ribbon_curves.first.size(); ++i)
			material_frames.push_back(ribbon_material_frame(i));
		
		return material_frames;
		
	};
	
	
	// ribbon material frames accessor
	Triad DiscreteRibbon::ribbon_material_frame(Size i) const {
		
		// related segments
		const Segment base_segment =
        CurveTopology::curve_segment(m_ribbon_curves.first, i);
		const Segment edge_segment =
        CurveTopology::curve_segment(m_ribbon_curves.second, i);
		
		// normal projection
		const Vector3 material_vector =
        Segment::normal_projection_axis(base_segment, edge_segment);
		
		// material frame
        Triad mat_frame;
        Vector3 v = base_segment.tangent().cross(material_vector);
        v.normalize();
        mat_frame.set_origin(base_segment.midpoint());
        mat_frame.set_axis(I, material_vector);
        mat_frame.set_axis(J, v);
        mat_frame.set_axis(K, base_segment.tangent());
		return mat_frame;

	};
		
	
}

