// LegoBender class
// Nicolas Clauvelin


// bending LEGO implementation


#ifndef emDNA_LegoBender_h
#define emDNA_LegoBender_h


#include <LegoProtein.h>


class LegoBender final : public LegoProtein {


public:

    // constructors
    LegoBender() = default;
    LegoBender(const LegoBender& lego_bender) = default;
    LegoBender(LegoBender&& lego_bender) = default;
    LegoBender(const Size& n_bp, const Real& bending_angle);

    // copy and move operators
    LegoBender& operator=(const LegoBender& lego_bender) = default;
    LegoBender& operator=(LegoBender&& lego_bender) = default;


private:

    // step parameters computation method
    void compute_bp_collection() override;

    // lego parameters
    Size m_n_bp;
    Real m_bending_angle = FLOAT_INIT;


};


#endif  // emDNA_LegoBender_h
