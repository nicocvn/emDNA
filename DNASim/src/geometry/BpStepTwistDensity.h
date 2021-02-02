// BpTwistDensity class
// Nicolas Clauvelin


// this class implements the calculation of the twist of "supercoiling" for
// a base-pair step as defined in:
// Two perspectives on the twist of DNA,
// L. A. Britton, W. K. Olson and I. Tobias
// Journal of Chemical Physics, 131 (24), 245101, 2009


#ifndef DNASim_BpStepTwistDensity_h
#define DNASim_BpStepTwistDensity_h


#include "DNASim_Includes.h"


namespace DNASim {


    class Segment;
    class Triad;
    typedef Triad BasePair;


	class BpTwistDensity {


	public:

        // static methods for twist density computation
        //
        // the step method returns three values:
        // - first is the twist density obtained with the short axis,
        // - second is obtained with the long axis,
        // - third is obtained with an averaged short/long axis
        static
        std::vector<Vector3>
        compute_twist_density(const std::vector<BasePair>& base_pairs);
        static
        Vector3
        compute_bp_step_twist_density(const BasePair& bp1,
                                      const BasePair& bp2,
                                      const Triad& vertex_frame1,
                                      const Triad& vertex_frame2);
        

    private:

        // vertex frame computation method
        static Triad compute_vertex_frame(const Segment& ingoing,
                                          const Segment& outgoing);

        // private constructors and copy operators
        // this class only exposes static methods
        BpTwistDensity() = delete;
        BpTwistDensity(const BpTwistDensity& bp_twist_density) = delete;
        BpTwistDensity(BpTwistDensity&& bp_twist_density) = delete;
        ~BpTwistDensity() = delete;

        // copy and move operators
        BpTwistDensity&
        operator=(const BpTwistDensity& bp_twist_density) = delete;
        BpTwistDensity&
        operator=(BpTwistDensity&& bp_twist_density) = delete;


    };


}


#endif  // DNASim_BpStepTwistDensity_h
