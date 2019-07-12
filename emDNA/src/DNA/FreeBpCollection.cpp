// FreeBpCollection class
// Nicolas Clauvelin


#include <BpCollectionGeometry.h>
#include <BpCollectionGradient.h>
#include <FreeBpCollection.h>


// class default constructor
FreeBpCollection::FreeBpCollection() : BpCollection_Interface() {};


// class constructor by copy
FreeBpCollection::FreeBpCollection(const FreeBpCollection& free_collection) :
BpCollection_Interface(free_collection) {};


// class destructor
FreeBpCollection::~FreeBpCollection() {};


// copy operator
FreeBpCollection& FreeBpCollection::operator=(const FreeBpCollection&
                                              free_collection) {
    BpCollection_Interface::operator=(free_collection);
    return *this;
};


// bp collection free dofs gradient
// this is the gradient used for the minimization
VectorN FreeBpCollection::bp_collection_free_dofs_gradient() const {

    return BpCollectionGradient::free_collection_gradient(m_bp_collection,
                                                          m_idx_mgr,
                                                          m_elec_flag);

};


// update method
void FreeBpCollection::update_with_new_inline_free_dofs(const VectorN&
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


// frozen steps rebuilding method
std::vector<BpStepDofs> FreeBpCollection::
rebuild_full_dofs_from_free_dofs(const std::vector<BpStepDofs>& free_dofs)
const {
    return rebuild_dofs_from_free_dofs(m_bp_collection.n_of_bp_steps(),
                                       free_dofs);
};
