// BpStepParams class
// Nicolas Clauvelin


// base-pair step parameters class

// base-pair step parameters units:
//  TILT    degree
//  ROLL    degree
//  TWIST   degree
//  SHIFT   angstroms
//  SLIDE   angstroms
//  RISE    angstroms


#ifndef emDNA_BpStepParams_h
#define emDNA_BpStepParams_h


#include <emDNA_Includes.h>


// typedef
using BpStepParameter = BpStepParameters;


class BpStepParams {


public:

    // constructors
    BpStepParams() = default;
    BpStepParams(const BpStepParams& bp_step_params) = default;
    BpStepParams(BpStepParams&& bp_step_params) = default;
    BpStepParams(const BasePair& bp1, const BasePair& bp2);
    BpStepParams(const VectorN& params_values);
    BpStepParams(const std::string& string_params);
    ~BpStepParams() = default;

    // copy and move operators
    BpStepParams& operator=(const BpStepParams& bp_step_params) = default;
    BpStepParams& operator=(BpStepParams&& bp_step_params) = default;

    // parameters accessors/modifiers
    inline const Real& value(const BpStepParameter& i) const {
        return m_params[i];
    };
    inline Real& value(const BpStepParameter& i) {
        return m_params[i];
    };

    // inline vector
    inline const VectorN& inline_vector() const {
        return m_params.inline_vector();
    };
    inline VectorN& inline_vector() {
        return m_params.inline_vector();
    };

    // base pair rebuild method
    BasePair rebuild_step_last_base_pair(const BasePair& bp) const;

    // static computation methods
    static
    BpStepParams bp_step_parameters(const BasePair& bp1, const BasePair& bp2);
    static
    std::vector<BasePair> rebuild_bps(const std::vector<BpStepParams>& prms,
                                      const BasePair& first_bp);
    static
    std::vector<BpStepParams> compute_params(const std::vector<BasePair> bps);

    
private:

    StepParameters m_params;


};


// const iterator typedef
typedef std::vector<BpStepParams>::const_iterator BpStepParamsConstIt;


#endif  // emDNA_BpStepParams_h
