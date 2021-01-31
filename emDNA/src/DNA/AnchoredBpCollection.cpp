// AnchoredBpCollection class
// Nicolas Clauvelin


#include "DNA/BpCollectionGeometry.h"
#include "DNA/BpCollectionGradient.h"
#include "DNA/AnchoredBpCollection.h"


// class default constructor
AnchoredBpCollection::AnchoredBpCollection() : BpCollection_Interface() {};


// class constructor by copy
AnchoredBpCollection::AnchoredBpCollection(const AnchoredBpCollection&
                                           anchored_collection) :
BpCollection_Interface(anchored_collection) {};


// class destructor
AnchoredBpCollection::~AnchoredBpCollection() {};


// copy operator
AnchoredBpCollection&
AnchoredBpCollection::operator=(const AnchoredBpCollection&
                                anchored_collection) {
    BpCollection_Interface::operator=(anchored_collection);
    return *this;
};


// bp collection free dofs accessors
// this is the set of dofs on which the minimization is performed
VectorN AnchoredBpCollection::bp_collection_inline_free_dofs() const {

    // inline free dofs without bc reduction
    const VectorN inline_free_dofs =
    BpCollection_Interface::bp_collection_inline_free_dofs();

    // the EED boundary conditions reduces "half" of the last free
    // step, that is, the last three free dofs are reduced
    const Size n_of_reduced_dofs =
    m_idx_mgr.n_of_free_dofs_variables()-
    BoundaryConditionsDofsTrimming::EEDCollection;

    // filling
    VectorN reduced_free_dofs(n_of_reduced_dofs, FLOAT_INIT);
    reduced_free_dofs.set_slice(0,
                                inline_free_dofs.slice(0, n_of_reduced_dofs));

    return reduced_free_dofs;

};


// bp collection free dofs scalings
// this is the scalings for the free dofs used in the minimization
VectorN AnchoredBpCollection::bp_collection_free_dofs_scalings() const {

    // free dofs scalings without bc reduction
    const VectorN free_dofs_scalings =
    BpCollection_Interface::bp_collection_free_dofs_scalings();

    // the EED boundary conditions reduces "half" of the last free
    // step, that is, the last three free dofs are reduced
    const Size n_of_reduced_dofs =
    m_idx_mgr.n_of_free_dofs_variables()-
    BoundaryConditionsDofsTrimming::EEDCollection;

    // filling
    VectorN reduced_dofs_scalings(n_of_reduced_dofs, FLOAT_INIT);
    reduced_dofs_scalings.set_slice(0,
                                    free_dofs_scalings.
                                    slice(0, n_of_reduced_dofs));

    return reduced_dofs_scalings;

};


// bp collection free dofs gradient
// this is the gradient used for the minimization
VectorN AnchoredBpCollection::bp_collection_free_dofs_gradient() const {

    return BpCollectionGradient::EED_collection_gradient(m_bp_collection,
                                                         m_idx_mgr,
                                                         m_elec_flag);

};


// update method
// there is no frozen steps in this collection
void AnchoredBpCollection::update_with_new_inline_free_dofs(const VectorN&
                                                            inline_free_dofs) {

    // dofs repacking
    std::vector<BpStepDofs> free_dofs =
    repack_inline_vector<BpStepDofs>(inline_free_dofs,
                                     emDNAConstants::StepParametersDim);

    // check if we have the expected number of packed step dofs
    // for imposed EED the number of free dofs should be the number of free
    // steps
    DS_ASSERT(free_dofs.size() == m_idx_mgr.n_of_free_steps(),
              "wrong size for the number of dofs");

    // rebuild the dofs and update the collection
    const std::vector<BpStepDofs> full_dofs =
    rebuild_full_dofs_from_free_dofs(free_dofs);
    m_bp_collection.update_collection_from_bp_step_dofs(full_dofs);

};


// frozen steps rebuilding method
std::vector<BpStepDofs> AnchoredBpCollection::
rebuild_full_dofs_from_free_dofs(const std::vector<BpStepDofs>&
                                 free_dofs) const {

    // bc reduced step
    const Size bc_step_index = m_idx_mgr.bc_reduced_step().collection_index();

    // first base pair
    const BasePair& first_bp = m_bp_collection.base_pair(0);

    // first step is to rebuild the dofs until the bc reduced step (excluded)
    const std::vector<BpStepDofs> dofs_until_bc_step =
    rebuild_dofs_from_free_dofs(bc_step_index, free_dofs);

    // now, we rebuild the base pairs until the bc step
    std::vector<BasePair> bps_until_bc_step =
    BpStepDofs::rebuild_bps(dofs_until_bc_step, first_bp);

    // the first bp is the last one of the previously computed set
    // the last bp is given by the base pair collection but we will only use
    // the translational dofs
    const BasePair& bc_step_first_bp = bps_until_bc_step[bc_step_index];
    const BasePair& bc_step_last_bp =
    m_bp_collection.base_pair(bc_step_index+1);

    // bc step dofs
    // this assumes that the last free dofs are corresponding to the bc step
    BpStepDofs bc_step_dofs(bc_step_first_bp, bc_step_last_bp);
    bc_step_dofs.value(TILTrad) = free_dofs.back().value(TILTrad);
    bc_step_dofs.value(ROLLrad) = free_dofs.back().value(ROLLrad);
    bc_step_dofs.value(TWISTrad) = free_dofs.back().value(TWISTrad);

    // we now have all the free dofs and can proceed with the full set of dofs
    std::vector<BpStepDofs> rebuild_free_dofs = free_dofs;
    rebuild_free_dofs.pop_back();
    rebuild_free_dofs.push_back(bc_step_dofs);

    return rebuild_dofs_from_free_dofs(m_bp_collection.n_of_bp_steps(),
                                       rebuild_free_dofs);

};
