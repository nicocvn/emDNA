// StepParameters class
// Nicolas Clauvelin


#include "maths/Vector3.h"
#include "maths/Matrix4.h"
#include "maths/Triad.h"
#include "dna/StepParameters.h"
#include "dna/BpGeometry.h"


// step parameters vector dimension
#define STEP_DIM 6


namespace DNASim {
	
	
	// class default constructor
	StepParameters::StepParameters() : m_parameters(STEP_DIM, FLOAT_INIT) {};


    // constructor from a initializer list
    StepParameters::StepParameters(const std::initializer_list<Real>& p) :
    StepParameters(VectorN(p)) {
        // size checking
		DS_ASSERT(m_parameters.size() == STEP_DIM,
				  "StepParameters(Vector<Real>) called with wrong size;"
				  " size=" << p.size());
    };
	

	// class constructor with initialization from a Vector
	StepParameters::StepParameters(const VectorN& v) :
        m_parameters(v) {
		// size checking
		DS_ASSERT(m_parameters.size() == STEP_DIM,
				  "StepParameters(Vector<Real>) called with wrong size;"
				  " size=" << v.size());
	};
	
	
	// class constructor with initialization from a string
	StepParameters::StepParameters(const std::string& s) : m_parameters(s) {
		// size checking
		DS_ASSERT(m_parameters.size() == STEP_DIM,
				  "StepParameters(std::string) called with wrong size;"
				  " size=" << VectorN(s).size());
	};
	
	
	// class constructor with initialization from bp triads
	StepParameters::StepParameters(const Triad& bp1, const Triad& bp2) :
	m_parameters() {
		*this = BpGeometry::bp_step_parameters_from_triads(bp1,bp2);
	};

		
	// step geometry methods
	Matrix4 StepParameters::step_matrix() const {
		return BpGeometry::step_matrix_from_parameters(*this);
	};
	Vector3 StepParameters::step_euler_zyz_angles() const {
		return BpGeometry::euler_zyz_angles_from_parameters(*this);
	};
	Vector3 StepParameters::step_local_displacement() const {
		return BpGeometry::local_step_displacement(*this);
	};
	Triad StepParameters::reconstruct_triad(const Triad& initial) const {
		return BpGeometry::reconstruct_bp_from_parameters(initial, *this);
	};


    // specialized print_vector function
    template <>
    void print_vector<StepParameters>(const std::vector<StepParameters>& data,
                                      const char& opening,
                                      const char& closing,
                                      const char& seprator) {
        Size n = data.size();
        std::cout << opening;
        for (Size i=0; i<n; ++i) {
            std::cout << data[i].inline_vector();
            if (i != n-1)
                std::cout << seprator;
        };
        std::cout << closing;
    };

		
}


#undef STEP_DIM
