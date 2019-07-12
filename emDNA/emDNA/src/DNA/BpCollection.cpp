// BpCollection class
// Nicolas Clauvelin


#include <BpCollection.h>


// instanciation from a set of base pairs
BpCollection
BpCollection::collection_from_base_pairs(const std::vector<BasePair>& bps) {

    // new collection
    BpCollection bp_collection;

    // update
    bp_collection.update_collection_from_base_pairs(bps);

    return bp_collection;

};


// instanciation from a set of bp step dofs
BpCollection
BpCollection::collection_from_bp_step_dofs(const std::vector<BpStepDofs>& dofs,
                                           const BasePair& first_bp) {

    // new collection
    BpCollection bp_collection;

    // update
    bp_collection.m_base_pairs.push_back(first_bp);
    bp_collection.update_collection_from_bp_step_dofs(dofs);

    return bp_collection;

};


// instanciation from a set of bp step params
BpCollection
BpCollection::
collection_from_bp_step_params(const std::vector<BpStepParams>& prms,
                               const BasePair& first_bp) {

    // new collection
    BpCollection bp_collection;

    // update
    bp_collection.m_base_pairs.push_back(first_bp);
    bp_collection.update_collection_from_bp_step_params(prms);

    return bp_collection;

};


// collection properties
Size BpCollection::n_of_base_pairs() const {
    return m_base_pairs.size();
};
Size BpCollection::n_of_bp_steps() const {
    return m_bp_step_dofs.size();
};


// single data acessors
const BasePair& BpCollection::base_pair(Size bp_index) const {
    return m_base_pairs[bp_index];
};
const BpStepParams& BpCollection::bp_step_params(Size step_index) const {
    return m_bp_step_params[step_index];
};
const BpStepDofs& BpCollection::bp_step_dofs(Size step_index) const {
    return m_bp_step_dofs[step_index];
};


// vector data accessors
const std::vector<BasePair>& BpCollection::base_pairs() const {
    return m_base_pairs;
};
const std::vector<BpStepParams>& BpCollection::bp_step_params() const {
    return m_bp_step_params;
};
const std::vector<BpStepDofs>& BpCollection::bp_step_dofs() const {
    return m_bp_step_dofs;
};


// inline data accessors
VectorN BpCollection::inline_all_bp_step_parameters() const {
    VectorN inline_params(m_bp_step_params.size()*
                          emDNAConstants::
                          StepParametersDim,
                          FLOAT_INIT);
    const Size n_steps = m_bp_step_params.size();
    for (Size i=0; i<n_steps; ++i) {
        inline_params.set_slice(emDNAConstants::StepParametersDim*i,
                                m_bp_step_params[i].inline_vector());
    };
    return inline_params;
};
VectorN BpCollection::inline_all_bp_step_dofs() const {
    VectorN inline_dofs(m_bp_step_dofs.size()*
                        emDNAConstants::
                        StepParametersDim,
                        FLOAT_INIT);
    const Size n_steps = m_bp_step_dofs.size();
    for (Size i=0; i<n_steps; ++i) {
        inline_dofs.set_slice(emDNAConstants::StepParametersDim*i,
                              m_bp_step_dofs[i].inline_vector());
    };
    return inline_dofs;
};


// frozen steps modifier/accessor
void BpCollection::set_frozen_steps_domains(const std::vector<SizePair>&
                                            frozen_parts) {
    m_frozen_steps_domains.clear();
    m_frozen_steps_domains = frozen_parts;
};
const std::vector<SizePair>& BpCollection::frozen_steps_domains() const {
    return m_frozen_steps_domains;
};


// sequence accesor/modifier
std::string BpCollection::collection_sequence() const {
    std::string full_sequence;
    auto end = m_step_sequences.end();
    for (auto it = m_step_sequences.begin(); it != end; ++it)
        full_sequence += Base::str(it->first_base());
    full_sequence += Base::str(m_step_sequences.back().last_base());
    return full_sequence;
};
void BpCollection::set_collection_sequence(const std::string& sequence) {
    m_step_sequences.clear();
    for (Size i=0; i<sequence.size()-1; ++i) {
        m_step_sequences.
        push_back(StepSequence(Base::
                               base_symbol_from_char(sequence[i]),
                               Base::
                               base_symbol_from_char(sequence[i+1])));
    };
};
void BpCollection::set_collection_dummy_sequence() {
    m_step_sequences.clear();
    for (Size i=0; i<n_of_bp_steps(); ++i)
        m_step_sequences.push_back(StepSequence(BaseSymbol::A,BaseSymbol::A));
};


// sequence-dependence model accessor/modifiers
std::string BpCollection::sequence_dependence_model() const {
    return std::string("step_model="+m_step_seqdep.model_name()
                       +
                       "+"
                       +"force_model="+m_fmat_seqdep.model_name());
};
void BpCollection::set_sequence_dependence_model(const std::string& model_name)
{
    m_step_seqdep.init_with_model(model_name);
    m_fmat_seqdep.init_with_model(model_name);
};
void BpCollection::
set_sequence_dependence_model(const StepParametersDB& steps_db,
                              const ForceConstantsDB& fmat_db) {
    m_step_seqdep = steps_db;
    m_fmat_seqdep = fmat_db;
};


// step sequence-dependent data accessors
const StepSequence& BpCollection::bp_step_sequence(Size step_index) const {
    return m_step_sequences[step_index];
};
const BpStepParams
BpCollection::bp_step_intrinsic_parameters(Size step_index) const {

    // StepParameters data
    const VectorN p =
    m_step_seqdep.intrinsic_bp_step_params(m_step_sequences[step_index]).
    inline_vector();

    // convert to BpStepParams
    return BpStepParams(p);

};
const MatrixN& BpCollection::bp_step_force_constants(Size step_index) const {
    return m_fmat_seqdep.force_constants(m_step_sequences[step_index]);
};


// update methods from base pairs
void BpCollection::
update_collection_from_base_pairs(const std::vector<BasePair>& bps) {

    // new base pairs
    m_base_pairs = bps;

    // orthogonalize all base pair frames
    for (BasePair& bp : m_base_pairs)
        bp.orthogonalize();

    // new step dofs
    m_bp_step_dofs = BpStepDofs::compute_dofs(bps);

    // new step params
    m_bp_step_params = BpStepParams::compute_params(bps);

};


// update methods from bp step parameters
void BpCollection::
update_collection_from_bp_step_params(const std::vector<BpStepParams>&
                                      bp_step_params) {

    // new step params
    m_bp_step_params = bp_step_params;

    // new base pairs
    m_base_pairs = BpStepParams::rebuild_bps(bp_step_params,
                                             m_base_pairs.front());

    // new step dofs
    m_bp_step_dofs = BpStepDofs::compute_dofs(m_base_pairs);

};


// update methods from bp step dofs
void BpCollection::
update_collection_from_bp_step_dofs(const std::vector<BpStepDofs>&
                                    bp_step_dofs) {

    // new step dofs
    m_bp_step_dofs = bp_step_dofs;

    // new base pairs
    m_base_pairs = BpStepDofs::rebuild_bps(bp_step_dofs, m_base_pairs.front());

    // new step params
    m_bp_step_params = BpStepParams::compute_params(m_base_pairs);

};

