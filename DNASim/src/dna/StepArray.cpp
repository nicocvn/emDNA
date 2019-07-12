// StepArray class
// Nicolas Clauvelin


#include <Vector3.h>
#include <VectorN.h>
#include <Triad.h>
#include <PRN_Integer.h>
#include <EnhancedString.h>
#include <Parser_x3DNA.h>
#include <StepArray.h>


namespace DNASim {
	

	// class constructor with size and default value initialization
	StepArray::StepArray(Size n, const StepParameters& p) :
	m_steps(n, p), m_generators() {
		
		// set generator matrices
		m_generators.reserve(n);
		for (Size i=0; i<n; ++i)
			m_generators.push_back(p.step_matrix());
		
	};

	
	// parameters accessor/modifier
	const StepParameters& StepArray::parameters(Size i) const {
		return m_steps[i];
	};
    const std::vector<StepParameters>& StepArray::all_parameters() const {
        return m_steps;
    };
	void StepArray::set_parameters(Size i, const StepParameters& p) {
		m_steps[i] = p;
		m_generators[i]	= p.step_matrix();
	};
	void StepArray::set_all_parameters(const StepParameters& p) {
		for (Size i=0, end=n_step(); i<end; ++i) {
			m_steps[i] = p;
			m_generators[i]	= p.step_matrix();
		};
	};
	
	
	// array management methods
	void StepArray::clear() {
		m_steps.clear();
		m_generators.clear();
	};
	void StepArray::append_step_parameters(const StepParameters& p) {
		m_steps.push_back(p);
		m_generators.push_back(p.step_matrix());
	};
	
	
	// single bp reconstruction method
	Triad StepArray::base_pair(Size i, const Triad& array_anchor) const {
		
		// first bp case
		if (i == 0)
			return array_anchor;
		
		// first bp frame matrix
		Matrix4 array_anchor_mat = array_anchor.matrix_representation();
		
		// accumulate rotation in the anchor bp matrix
        for (Size idx=0; idx<i; ++idx)
            array_anchor_mat *= m_generators[idx];
		
		return Triad(array_anchor_mat);

	};
	
	
	// set of bps reconstruction method
	std::vector<Triad> StepArray::base_pairs(Size first, Size last,
											 const Triad& array_anchor) const {

        // data container
		std::vector<Triad> bps;
		bps.reserve(last-first+1);
		
		// reconstruct first bp of the set
		Triad first_bp = base_pair(first, array_anchor);
		bps.push_back(first_bp);
		
		// reconstruct following bps
		Matrix4 mat = first_bp.matrix_representation();
		std::vector<Matrix4>::const_iterator it, end;
        for (Size i=first; i<last; ++i) {
			mat *= m_generators[i];
			bps.push_back(Triad(mat));
		};
		
		return bps;
		
	};
	
	
	// all bps reconstruction method
	std::vector<Triad> StepArray::all_base_pairs(const Triad& array_anchor)
	const {
        return base_pairs(0, n_bp()-1, array_anchor);
	};


    // static method for x3DNA format output
    std::string StepArray::x3DNA_form_output() const {

        // x3DNA data with dummy sequence
        Data_x3DNA<StepParameters> data;
        data._data = m_steps;
        char* seq = (char*)std::malloc(m_steps.size()+1);
        std::memset(seq, 'A', m_steps.size()+1);
        data._sequence = std::string(seq);
        std::free(seq);
        seq = NULL;

        // x3DNA parser
        std::stringstream ss;
        Parser_x3DNA::create_step_parameters_output(ss, data);
        return ss.str();

    };


};
