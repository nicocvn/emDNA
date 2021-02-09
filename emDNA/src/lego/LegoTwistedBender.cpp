// LegoTwistedBender class
// Nicolas Clauvelin


#include "lego/LegoTwistedBender.h"


// class constructor with full initialization
LegoTwistedBender::LegoTwistedBender(const Size& n_bp,
                                     const Real& added_twist,
                                     const Real& bending_angle) :
LegoProtein(),
m_n_bp(n_bp), m_added_twist(added_twist), m_bending_angle(bending_angle) {
    compute_bp_collection();
};


// step parameters computation method
void LegoTwistedBender::compute_bp_collection() {

    // IdealDNA force field
    // this force field is uniform so we can use dummy sequence
    StepParametersDB model("IdealDNA");

    // default intrinsic step parameters
    BpStepParams p0(model.
                    intrinsic_bp_step_params(StepSequence(BaseSymbol::A,
                                                          BaseSymbol::A)).
                    inline_vector());
    const Real& tw0 = p0.value(TWIST)*DEG_2_RAD;

    // sizing
    const Size n_steps = m_n_bp-1;

    // filling
    std::vector<BpStepParams> prms;
    prms.reserve(n_steps);
    for (Size i=0; i<n_steps; ++i) {
        VectorN pvec = {
            m_bending_angle*std::sin(i*(tw0+m_added_twist*DEG_2_RAD)),
            m_bending_angle*std::cos(i*(tw0+m_added_twist*DEG_2_RAD)),
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


