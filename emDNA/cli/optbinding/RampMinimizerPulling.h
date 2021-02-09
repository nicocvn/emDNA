// RampMinimizerPulling class
// Nicolas Clauvelin


// pulling binding ramp minimizer


#ifndef emDNA_RampMinimizerPulling_h
#define emDNA_RampMinimizerPulling_h


#include <RampMinimizer.h>


class RampMinimizerPulling final : public RampMinimizer {


public:

    // constructors
    RampMinimizerPulling() = default;
    RampMinimizerPulling(const Vector3& pulling_force);
    RampMinimizerPulling(const RampMinimizerPulling& ramp_min) = default;
    RampMinimizerPulling(RampMinimizerPulling&& ramp_min) = default;
    ~RampMinimizerPulling() = default;

    // copy and move operators
    RampMinimizerPulling&
    operator=(const RampMinimizerPulling& ramp_min) = default;
    RampMinimizerPulling& operator=(RampMinimizerPulling&& ramp_min) = default;

    // minimization method
    bool binding_ramp_minimization() override;

    // pulling force accessor
    const Vector3& pulling_force() const;


private:

    std::vector<std::string> headers_list() const override {
        return {
            "ENERGY", "|GRAD_NORM|", "ITERATIONS", "RETURN_CODE",
            "CONTRIB_1", "VAL_1", "CONTRIB_2", "VAL_2"
        };
    };

    Vector3 m_pulling_force;
    

};


#endif  // emDNA_RampMinimizerPulling_h
