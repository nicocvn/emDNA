// StepArray class
// Nicolas Clauvelin


// base pair step container implementation


#ifndef DNASim_StepArray_h
#define DNASim_StepArray_h


#include "dna/StepParameters.h"
#include "maths/Matrix4.h"


namespace DNASim {
	
	
	class Triad;
	
	
	class StepArray {
		
		
	public:
		
		// constructors
		StepArray() = default;
		StepArray(Size n, const StepParameters& p);
		StepArray(const StepArray& steps_array) = default;
        StepArray(StepArray&& steps_array) = default;
		~StepArray() = default;
		
		// copy and move operators
		StepArray& operator=(const StepArray& steps_array) = default;
        StepArray& operator=(StepArray&& steps_array) = default;
		
		// array properties
		inline Size	n_step() const { return m_steps.size(); };
		inline Size	n_bp() const { return m_steps.size()+1; };
		
		// parameters accessor/modifier
		const StepParameters& parameters(Size i) const;
        const std::vector<StepParameters>& all_parameters() const;
		void set_parameters(Size i, const StepParameters& p);
		void set_all_parameters(const StepParameters& p);
		
		// array management methods
		void clear();
		void append_step_parameters(const StepParameters& p);
				
		// bp reconstruction methods
        // in base_pairs(first,last,anchor) last is inclusive
		Triad base_pair(Size i, const Triad& array_anchor) const;
		std::vector<Triad> base_pairs(Size first, Size last,
									  const Triad& array_anchor) const;
		std::vector<Triad> all_base_pairs(const Triad& array_anchor) const;
		
		// array iterators
		inline std::vector<StepParameters>::const_iterator begin() const {
			return m_steps.begin();
		};
		inline std::vector<StepParameters>::const_iterator end() const {
			return m_steps.end();
		};

        // static method for x3DNA format output
        std::string x3DNA_form_output() const;
		
		
	private:
		
		std::vector<StepParameters>	m_steps;
		std::vector<Matrix4> m_generators;
		
		
	};
	
	
}


#endif
