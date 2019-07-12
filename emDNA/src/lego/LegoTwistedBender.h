// LegoTwistedBender class
// Nicolas Clauvelin


#ifndef emDNA_LegoTwistedBender_h
#define emDNA_LegoTwistedBender_h


#include <LegoProtein.h>


class LegoTwistedBender final : public LegoProtein {


public:

    // constructors
    LegoTwistedBender() = default;
    LegoTwistedBender(const LegoTwistedBender& lego_twister) = default;
    LegoTwistedBender(LegoTwistedBender&& lego_twister) = default;
    LegoTwistedBender(const Size& n_bp,
                      const Real& added_twist,
                      const Real& bending_angle);

    // copy and move operators
    LegoTwistedBender& operator=(const LegoTwistedBender&
                                 lego_twister) = default;
    LegoTwistedBender& operator=(LegoTwistedBender&& lego_twister) = default;


private:

    // step parameters computation method
    void compute_bp_collection() override;

    // lego parameters
    Size m_n_bp;
    Real m_bending_angle = FLOAT_INIT;
    Real m_added_twist = FLOAT_INIT;

};


#endif  // emDNA_LegoTwistedBender_h
