// DiscreteRibbon class
// Nicolas Clauvelin


// implementation of methods and formula described in:
// Characterization of the Geometry and Topology of DNA Pictured As a Discrete
// Collection of Atoms
// N. Clauvelin, W. K. Olson, I. Tobias
// JCTC, 2012, DOI: 10.1021/ct200657e


#ifndef DNASim_DiscreteRibbon_h
#define DNASim_DiscreteRibbon_h


#include <DNASim_Includes.h>


namespace DNASim {
	
	
	class Triad;
	typedef Vector3 Vertex;
	typedef std::vector<Vertex> CurvePoints;
	
	
	class DiscreteRibbon {
		
		
	public:
		
		// constructors
		DiscreteRibbon(const CurvePoints& base_curve,
					   const CurvePoints& edge_reference_curve);
		~DiscreteRibbon();
		
		// ribbon curves accessors
		const CurvePoints&	ribbon_base_curve() const;
		const CurvePoints&	ribbon_edge_reference_curve() const;
		
		// ribbon topology methods
		Real ribbon_writhing_number() const;
		Real ribbon_linking_number() const;
		Real ribbon_total_twist() const;
		
		// writhe gradient method
		std::vector<Vector3> ribbon_writhe_gradient() const;
		
		// ribbon vertex property
		Real vertex_curvature(Size i) const;
		Real vertex_twist(Size i) const;
		Real vertex_length(Size i) const;
		
		// ribbon material frames accessors
		std::vector<Triad>	ribbon_material_frames() const;
		Triad				ribbon_material_frame(Size i) const;
		
		
	private:
		
		std::pair<CurvePoints,CurvePoints> m_ribbon_curves;
		
		
	};
	
	
}


#endif	// DNASim_DiscreteRibbon_h
