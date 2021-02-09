// RampMinimizerEEDR class
// Nicolas Clauvelin


// EEDR binding ramp minimizer


#ifndef emDNA_RampMinimizerEEDR_h
#define emDNA_RampMinimizerEEDR_h


#include <RampMinimizer.h>


class RampMinimizerEEDR final : public RampMinimizer {


public:

    // constructors
    RampMinimizerEEDR() = default;
    RampMinimizerEEDR(const RampMinimizerEEDR& ramp_min) = default;
    RampMinimizerEEDR(RampMinimizerEEDR&& ramp_min) = default;
    ~RampMinimizerEEDR() = default;

    // copy and move operators
    RampMinimizerEEDR& operator=(const RampMinimizerEEDR& ramp_min) = default;
    RampMinimizerEEDR& operator=(RampMinimizerEEDR&& ramp_min) = default;

    // minimization method
    bool binding_ramp_minimization() override;
    

};


#endif  // emDNA_RampMinimizerEEDR_h
