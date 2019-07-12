// LegoTwister class
// Nicolas Clauvelin


// twisting LEGO implementation


#ifndef emDNA_LegoTwister_h
#define emDNA_LegoTwister_h


#include <LegoProtein.h>


class LegoTwister final : public LegoProtein {


public:

    // constructors
    LegoTwister() = default;
    LegoTwister(const LegoTwister& lego_twister) = default;
    LegoTwister(LegoTwister&& lego_twister) = default;
    LegoTwister(const Size& n_bp, const Real& added_twist);

    // copy and move operators
    LegoTwister& operator=(const LegoTwister& lego_twister) = default;
    LegoTwister& operator=(LegoTwister&& lego_twister) = default;


private:

    // step parameters computation method
    void compute_bp_collection() override;

    // lego parameters
    Size m_n_bp;
    Real m_added_twist = FLOAT_INIT;
    
    
};




#endif  // emDNA_LegoTwister_h
