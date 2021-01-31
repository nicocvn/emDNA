// BpCollection_Interface class
// Nicolas Clauvelin


#include "DNA/BpCollectionElasticEnergy.h"
#include "DNA/BpCollectionElectrostaticEnergy.h"
#include "DNA/BpCollectionGradient.h"
#include "DNA/BpCollection_Interface.h"


// bp collection accessor/modifier
const BpCollection& BpCollection_Interface::bp_collection() const {
    return m_bp_collection;
};
void BpCollection_Interface::set_bp_collection(const BpCollection&
                                               bp_collection) {

    // new bp collection
    m_bp_collection = bp_collection;

    // index manager setup
    m_idx_mgr.set_n_of_bp_step(bp_collection.n_of_bp_steps());
    m_idx_mgr.set_frozen_steps_indexes(bp_collection.frozen_steps_domains());

};


// electrostatic interactions flag modifiers
void BpCollection_Interface::toggle_on_electrostatics() {
    m_elec_flag = true;
};
void BpCollection_Interface::toggle_off_electrostatics() {
    m_elec_flag = false;
};


// bp collection energy
// this is the objective function
Real BpCollection_Interface::bp_collection_energy() const {

    Real energy = FLOAT_INIT;

    // elastic energy computed over all free steps
    energy +=
    BpCollectionElasticEnergy::elastic_energy(m_bp_collection, m_idx_mgr);

    // electrostatics energy
    if (m_elec_flag)
        energy +=
        BpCollectionElectrostaticEnergy::electrostatic_energy(m_bp_collection);

    return energy;

};


// bp collection energy split
// returns the list of contributions to the energy
std::vector<EnergyContrib> BpCollection_Interface::energy_contribs() const {

    std::vector<EnergyContrib> contribs;

    // elastic energy contrib is always here by default
    contribs.
    push_back(EnergyContrib("elastic",
                            BpCollectionElasticEnergy::
                            elastic_energy(m_bp_collection, m_idx_mgr)));

    // electrostatic energy contrib only if accounted for
    if (m_elec_flag)
        contribs.
        push_back(EnergyContrib("electrostatic",
                                BpCollectionElectrostaticEnergy::
                                electrostatic_energy(m_bp_collection)));

    // energy contribs
    return contribs;

};


// bp collection free dofs accessor
// this is the set of dofs on which the minimization is performed
VectorN BpCollection_Interface::bp_collection_inline_free_dofs() const {

    // inline free dofs indices
    const std::vector<Size> free_dofs_indices =
    m_idx_mgr.inline_free_dofs_indices();

    // all step dofs
    const VectorN all_dofs = m_bp_collection.inline_all_bp_step_dofs();

    // selection
    VectorN inline_free_dofs(free_dofs_indices.size(), FLOAT_INIT);
    for (Size i=0, end = free_dofs_indices.size(); i<end; ++i)
        inline_free_dofs[i] = all_dofs[free_dofs_indices[i]];

    return inline_free_dofs;

};


// bp collection free dofs scalings
// this is the scalings for the free dofs used in the minimization
VectorN BpCollection_Interface::bp_collection_free_dofs_scalings() const {

    // container
    VectorN all_dofs_scalings(m_idx_mgr.n_of_dofs_variables(), FLOAT_INIT);

    // loop over all steps
    const Size n_steps = m_bp_collection.n_of_bp_steps();
    for (Size i=0; i<n_steps; ++i) {

        // put the step dofs in the global container
        VectorN dofs_scalings(6, Real(1));
        dofs_scalings[0] = DEG_2_RAD;
        dofs_scalings[1] = DEG_2_RAD;
        dofs_scalings[2] = DEG_2_RAD;
        const Size global_index = m_idx_mgr.dof_global_coordinate(i, TILTrad);
        all_dofs_scalings.set_slice(global_index, dofs_scalings);

    };

    // inline free steps collection indexes
    const std::vector<Size> free_indexes = m_idx_mgr.inline_free_dofs_indices();

    // selection over the free steps
    VectorN free_dofs_scalings(free_indexes.size(), FLOAT_INIT);
    for (Size end = free_indexes.size(), i=0; i<end; ++i)
        free_dofs_scalings[i] = all_dofs_scalings[free_indexes[i]];

    return free_dofs_scalings;

};


// bp collection elastic energy gradient
// this is the gradient of the elastic energy only
VectorN BpCollection_Interface::bp_collection_elastic_energy_gradient() const {
    return BpCollectionGradient::free_collection_gradient(m_bp_collection,
                                                          m_idx_mgr,
                                                          m_elec_flag);
};


// class default constructor
BpCollection_Interface::BpCollection_Interface() :
m_bp_collection(),
m_idx_mgr(),
m_elec_flag(false) {};


// class constructor by copy
BpCollection_Interface::BpCollection_Interface(const BpCollection_Interface&
                                               bp_collection_interface) :
m_bp_collection(bp_collection_interface.m_bp_collection),
m_idx_mgr(bp_collection_interface.m_idx_mgr),
m_elec_flag(bp_collection_interface.m_elec_flag) {};


// copy operator
BpCollection_Interface&
BpCollection_Interface::operator=(const BpCollection_Interface&
                                  bp_collection_interface) {

    // set bp collection (index manager is setup automatically)
    set_bp_collection(bp_collection_interface.m_bp_collection);

    // elec flag
    m_elec_flag = bp_collection_interface.m_elec_flag;

    return *this;
    
};


// method for rebuilding dofs from free dofs
// dofs for the step in the range [start, end) are rebuild
std::vector<BpStepDofs> BpCollection_Interface::
rebuild_dofs_from_free_dofs(const Size& range_end,
                            const std::vector<BpStepDofs>& free_dofs) const {

    // dofs container
    std::vector<BpStepDofs> dofs(range_end);

    // free steps iterators
    const auto free_begin = m_idx_mgr.free_step_indexes_begin();
    const auto free_end = m_idx_mgr.free_step_indexes_end();

    // fill in the known dofs from the free steps
    std::for_each(free_begin, free_end,
                  [&dofs, &free_dofs, &range_end](const IndexPair& idx){

                      // check if the free step is in the range
                      if (idx.collection_index() < range_end)
                          dofs[idx.collection_index()] =
                          free_dofs[idx.local_index()];

                  });

    // range first base pair (anchoring bp)
    const BasePair& first_bp = m_bp_collection.base_pair(0);

    // loop over the frozen steps
    const auto frozen_begin = m_idx_mgr.frozen_step_indexes_begin();
    const auto frozen_end = m_idx_mgr.frozen_step_indexes_end();
    for (auto frozen_it = frozen_begin; frozen_it != frozen_end; ++frozen_it) {

        // collection index for the frozen step
        const Size& coll_idx = frozen_it->collection_index();

        // if the frozen step is outside of the range we can break
        if (coll_idx >= range_end)
            break;

        // extract dofs until the frozen step
        const std::vector<BpStepDofs> tmp_dofs(dofs.begin(),
                                               dofs.begin()+coll_idx);

        // rebuild all base pairs until the frozen step
        std::vector<BasePair> bps = BpStepDofs::rebuild_bps(tmp_dofs,
                                                                  first_bp);

        // select first base pair of the frozen step
        const BasePair& frozen_first_bp = bps.back();

        // frozen step parameters
        // these parameters cannot be altered
        const StepParameters frozen_p(m_bp_collection.
                                      bp_step_params(coll_idx).inline_vector());

        // compute last base pair of the frozen step
        BasePair frozen_last_bp =
        frozen_p.reconstruct_triad(frozen_first_bp);

        // compute the dofs corresponding to the frozen steps
        dofs[coll_idx] = BpStepDofs::bp_step_dofs(frozen_first_bp,
                                                  frozen_last_bp);
        
    };

    return dofs;

};

