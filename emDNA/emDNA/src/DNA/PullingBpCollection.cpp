// PulingBpCollection class
// Nicolas Clauvelin


#include <BpCollectionGeometry.h>
#include <BpCollectionGradient.h>
#include <PullingBpCollection.h>


// class default constructor
PullingBpCollection::PullingBpCollection() :
BpCollection_Interface(),
m_pulling_force(),
m_rest_state_endpoint() {};


// class constructor by copy
PullingBpCollection::PullingBpCollection(const PullingBpCollection&
                                         pulled_collection) :
BpCollection_Interface(pulled_collection),
m_pulling_force(pulled_collection.m_pulling_force),
m_rest_state_endpoint(pulled_collection.m_rest_state_endpoint) {};


// class destructor
PullingBpCollection::~PullingBpCollection() {};


// copy operator
PullingBpCollection& PullingBpCollection::operator=(const PullingBpCollection&
                                                    pulled_collection) {
    BpCollection_Interface::operator=(pulled_collection);
    m_pulling_force = pulled_collection.m_pulling_force;
    m_rest_state_endpoint = pulled_collection.m_rest_state_endpoint;
    return *this;
};


// bp collection modifier
void PullingBpCollection::set_bp_collection(const BpCollection& bp_collection) {

    // base class modifier
    BpCollection_Interface::set_bp_collection(bp_collection);

    // rest state endpoint
    m_rest_state_endpoint =
    bp_collection.base_pair(bp_collection.n_of_base_pairs()-1).origin();

};


// bp collection energy
// this is the objective function
// the energy is always calculated on the set of free steps
Real PullingBpCollection::bp_collection_energy() const {

    Real E = BpCollection_Interface::bp_collection_energy();

    // pulling term
    E += pulling_energy();

    return E;

};


// bp collection energy contribs
// returns the list of contributions to the energy
std::vector<EnergyContrib> PullingBpCollection::energy_contribs() const {

    return {
        EnergyContrib("elastic",
                      BpCollection_Interface::bp_collection_energy()),
        EnergyContrib("pulling",
                      pulling_energy())
    };

};


// bp collection free dofs gradient
// this is the gradient used for the minimization
VectorN PullingBpCollection::bp_collection_free_dofs_gradient() const {

    return
    BpCollectionGradient::
    pulling_collection_gradient(m_bp_collection,
                                m_idx_mgr,
                                m_pulling_force,
                                m_elec_flag);

};


// update method
void PullingBpCollection::update_with_new_inline_free_dofs(const VectorN&
                                                           inline_free_dofs) {

    // dofs repacking
    const std::vector<BpStepDofs> free_dofs =
    repack_inline_vector<BpStepDofs>(inline_free_dofs,
                                     emDNAConstants::StepParametersDim);

    // check if we have the expected number of packed step dofs
    DS_ASSERT(free_dofs.size() == m_idx_mgr.n_of_free_steps(),
              "wrong size for the number of dofs");

    // if there is no frozen steps we can update the collection
    if (m_idx_mgr.n_of_frozen_steps() == 0) {
        m_bp_collection.update_collection_from_bp_step_dofs(free_dofs);
        return;
    };

    // if there are frozen steps we need to rebuild the set of all dofs
    const std::vector<BpStepDofs> full_dofs =
    rebuild_full_dofs_from_free_dofs(free_dofs);
    m_bp_collection.update_collection_from_bp_step_dofs(full_dofs);

};


// pulling force accessor/modifier
const Vector3& PullingBpCollection::pulling_force() const {
    return m_pulling_force;
};
void PullingBpCollection::set_pulling_force(const Vector3& f) {
    m_pulling_force = f;
};


// pulling energy method
Real PullingBpCollection::pulling_energy() const {

    // pulling contrib
    const Vector3& endpoint =
    m_bp_collection.base_pair(m_bp_collection.n_of_base_pairs()-1).origin()
    -
    m_rest_state_endpoint;

    return -m_pulling_force.dot(endpoint)*emDNAConstants::PullingForceScaling;

};


// frozen steps rebuilding method
std::vector<BpStepDofs> PullingBpCollection::
rebuild_full_dofs_from_free_dofs(const std::vector<BpStepDofs>& free_dofs)
const {
    return rebuild_dofs_from_free_dofs(m_bp_collection.n_of_bp_steps(),
                                       free_dofs);
};
