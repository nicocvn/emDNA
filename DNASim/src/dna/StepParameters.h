// StepParameters class
// Nicolas Clauvelin


// step parameters implementation (built on on VectorN class)


#ifndef DNASim_StepParameters_h
#define DNASim_StepParameters_h


#include "maths/VectorN.h"


// enum type for step parameters
enum BpStepParameters {TILT=0, ROLL=1, TWIST=2,SHIFT=3, SLIDE=4, RISE=5};


namespace DNASim {
    

	class Vector3;
	class Matrix4;
	class Triad;
	
	
	class StepParameters {
		
		
	public:
		
		// constructors
		// default constructor creates a zero step parameters vector
		StepParameters();
		StepParameters(const StepParameters& p) = default;
        StepParameters(StepParameters&& p) = default;
        StepParameters(const std::initializer_list<Real>& p);
        explicit StepParameters(const VectorN& v);
		explicit StepParameters(const std::string& s);
		explicit StepParameters(const Triad& bp1, const Triad& bp2);
		~StepParameters() = default;
		
		// copy and move operators
		StepParameters& operator=(const StepParameters& p) = default;
        StepParameters& operator=(StepParameters&& p) = default;
		
		// inline vector accessor/modifier
        inline const VectorN& inline_vector() const {
            return m_parameters;
        };
        inline VectorN& inline_vector() {
            return m_parameters;
        };
		
		// parameters accessors/modifiers
		inline const Real& operator[](enum BpStepParameters p) const {
            return m_parameters[p];
        };
		inline Real& operator[](enum BpStepParameters p) {
            return m_parameters[p];
        };
		
		// step geometry methods
		// the local displacement is expressed in a canonical triad (needs to 
		// be expressed with respect to the corresponding base pair)
		Matrix4	step_matrix() const;
		Vector3	step_euler_zyz_angles() const;
		Vector3	step_local_displacement() const;
		Triad reconstruct_triad(const Triad& initial) const;
		
		
	private:
		
		VectorN m_parameters;
		

	};
    

    // specialized print_vector function
    template <>
    void print_vector<StepParameters>(const std::vector<StepParameters>& data,
                                      const char& opening,
                                      const char& closing,
                                      const char& seprator);


}


#endif  // DNASim_StepParameters_h
