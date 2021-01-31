// PRN_Real class
// Nicolas Clauvelin


// static class for random floats generation
// generate a (real) number in the range [low;up)


#ifndef DNASim_PRN_Real_h
#define DNASim_PRN_Real_h


#include "DNASim_Includes.h"


namespace DNASim {
	

    template <RandomSeed flag>
	class PRN_Real {
		
		
	public:
		
		// generate in [0;1) method
		static Real generate() {
            std::uniform_real_distribution<Real> distrib(Real(0), Real(1));
			return distrib(engine());
		};
		
		// generate with bounds method
        // returns a number in the range [low,up)
		static Real generate_with_bounds(Real low, Real up) {
            std::uniform_real_distribution<Real> distrib(low, up);
			return distrib(engine());
		};
		
		
	private:
		
		PRN_Real() = delete;
        ~PRN_Real() = delete;
		PRN_Real(const PRN_Real& prn_real) = delete;
        PRN_Real(PRN_Real&& prn_real) = delete;
        PRN_Real& operator=(const PRN_Real& prn_real) = delete;
        PRN_Real& operator=(PRN_Real&& prn_real) = delete;

        // static engine
        static std::mt19937& engine() {
            static std::random_device rdev;
            if (flag == RandomSeed::Yes) {
                static std::mt19937 m_engine(rdev());
                return m_engine;
            };
            static std::mt19937 m_engine;
            return m_engine;
        };
        

	};
	
	
}


#endif	// DNASim_PRN_Real_h
