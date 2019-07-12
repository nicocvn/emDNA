// PRN_Integer class
// Nicolas Clauvelin


// static class for random integers generation
// returns an integer in the range [low;up]


#ifndef DNASim_PRN_Integer_h
#define DNASim_PRN_Integer_h


#include <DNASim_Includes.h>


namespace DNASim {
	

    template <RandomSeed flag>
	class PRN_Integer {
		
		
	public:
		
		// generate with bounds method
        // returns a number in the range [low,up]
		static Integer generate_with_bounds(Integer low, Integer up) {
            DS_ASSERT((1+up-low) > Integer(0),
                      "bad bounds for number generation");
            std::uniform_int_distribution<Integer> distrib(low,up);
            return distrib(engine());
		};


	private:
		
		PRN_Integer() = delete;
        ~PRN_Integer() = delete;
		PRN_Integer(const PRN_Integer& prn_int) = delete;
        PRN_Integer(PRN_Integer&& prn_int) = delete;
		PRN_Integer& operator=(const PRN_Integer& prn_int) = delete;
        PRN_Integer& operator=(PRN_Integer&& prn_int) = delete;

        // static engine
        static std::mt19937& engine() {
            static std::random_device m_rdev;
            if (flag == RandomSeed::Yes) {
                static std::mt19937 m_engine(m_rdev());
                return m_engine;
            };
            static std::mt19937 m_engine;
            return m_engine;
        };
        
		
	};
	
}


#endif	// DNASim_PRN_Integer_h
