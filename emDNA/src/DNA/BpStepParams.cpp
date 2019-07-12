// BpStepParams class
// Nicolas Clauvelin


#include <BpStepParams.h>


// class constructor with initialization from two base pairs
BpStepParams::BpStepParams(const BasePair& bp1, const BasePair& bp2) :
m_params(bp1, bp2) {};



// class constructor with initialization from a vector
BpStepParams::BpStepParams(const VectorN& params_values) :
m_params(params_values) {};


// class constructor with initialization from a string
BpStepParams::BpStepParams(const std::string& string_params) :
m_params(string_params) {};


// base pair rebuild method
BasePair BpStepParams::rebuild_step_last_base_pair(const BasePair& bp) const {
    BasePair last_bp = m_params.reconstruct_triad(bp);
    last_bp.orthogonalize();
    return last_bp;
};


// static computation methods
BpStepParams BpStepParams::bp_step_parameters(const BasePair& bp1,
                                              const BasePair& bp2) {

    BpStepParams p;
    p.m_params = StepParameters(bp1, bp2);

    return p;

};
std::vector<BasePair>
BpStepParams::rebuild_bps(const std::vector<BpStepParams>& prms,
                          const BasePair& first_bp) {

    // number of base pairs
    const Size n_bp = prms.size() + 1;

    // allocation and first base pair
    std::vector<BasePair> base_pairs;
    base_pairs.reserve(n_bp);
    base_pairs.push_back(first_bp);

    // rebuild
    BpStepParamsConstIt step_end = prms.end();
    for (auto step_it = prms.begin(); step_it != step_end; ++step_it) {
        base_pairs.push_back(step_it->
                             rebuild_step_last_base_pair(base_pairs.back()));
    };

    return base_pairs;

};
std::vector<BpStepParams>
BpStepParams::compute_params(const std::vector<BasePair> bps) {

    // number of steps
    const Size n_step = bps.size()-1;

    // allocation
    std::vector<BpStepParams> params;
    params.reserve(n_step);

    // computation
    BasePairConstIt bp_end = bps.end();
    for (auto bp_it = bps.begin(); bp_it != bp_end-1; ++bp_it)
        params.push_back(BpStepParams::bp_step_parameters(*(bp_it),
                                                          *(bp_it+1)));

    return params;

};
