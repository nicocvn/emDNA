// LegoTwister class
// Nicolas Clauvelin


#include "lego/LegoTwister.h"


// class constructor with full initialization
LegoTwister::LegoTwister(const Size& n_bp, const Real& added_twist) :
LegoProtein(), m_n_bp(n_bp), m_added_twist(added_twist) {
    compute_bp_collection();
};


// step parameters computation method
void LegoTwister::compute_bp_collection() {

    // IdealDNA force field
    // this force field is uniform so we can use dummy sequence
    StepParametersDB model("IdealDNA");

    // default intrinsic step parameters
    BpStepParams p0(model.
                    intrinsic_bp_step_params(StepSequence(BaseSymbol::A,
                                                          BaseSymbol::A)).
                    inline_vector());

    // sizing
    const Size n_steps = m_n_bp-1;

    // filling
    std::vector<BpStepParams> prms;
    prms.reserve(n_steps);
    for (Size i=0; i<n_steps; ++i) {
        VectorN pvec = {
            FLOAT_INIT,
            FLOAT_INIT,
            m_added_twist,
            FLOAT_INIT,
            FLOAT_INIT,
            FLOAT_INIT
        };
        pvec += p0.inline_vector();
        prms.push_back(pvec);
    };

    // collection
    m_lego_coll = BpCollection::collection_from_bp_step_params(prms,
                                                               BasePair());

};

