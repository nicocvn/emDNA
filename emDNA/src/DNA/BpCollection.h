// BpCollection class
// Nicolas Clauvelin


// base pair collection container implementation


#ifndef emDNA_BpCollection_h
#define emDNA_BpCollection_h


#include "DNA/BpStepParams.h"
#include "DNA/BpStepDofs.h"


class BpCollection {


public:

    // instantiation methods
    static BpCollection
    collection_from_base_pairs(const std::vector<BasePair>& bps);
    static BpCollection
    collection_from_bp_step_dofs(const std::vector<BpStepDofs>& dofs,
                                 const BasePair& first_bp);
    static BpCollection
    collection_from_bp_step_params(const std::vector<BpStepParams>& prms,
                                   const BasePair& first_bp);

    // constructors
    BpCollection() = default;
    BpCollection(const BpCollection& bp_collection) = default;
    BpCollection(BpCollection&& bp_collection) = default;
    ~BpCollection() = default;

    // copy and move operators
    BpCollection& operator=(const BpCollection& bp_collection) = default;
    BpCollection& operator=(BpCollection&& bp_collection) = default;

    // collection properties
    Size n_of_base_pairs() const;
    Size n_of_bp_steps() const;

    // single data acessors
    const BasePair& base_pair(Size bp_index) const;
    const BpStepParams& bp_step_params(Size step_index) const;
    const BpStepDofs& bp_step_dofs(Size step_index) const;

    // vector data accessors
    const std::vector<BasePair>& base_pairs() const;
    const std::vector<BpStepParams>& bp_step_params() const;
    const std::vector<BpStepDofs>& bp_step_dofs() const;

    // inline data accessors
    VectorN inline_all_bp_step_parameters() const;
    VectorN inline_all_bp_step_dofs() const;

    // frozen steps modifier/accessor
    // the domains are inclusive [lower, upper]
    void set_frozen_steps_domains(const std::vector<SizePair>& frozen_parts);
    const std::vector<SizePair>& frozen_steps_domains() const;

    // sequence accesor/modifier
    std::string collection_sequence() const;
    void set_collection_sequence(const std::string& sequence);
    void set_collection_dummy_sequence();

    // sequence-dependence model accessor/modifiers
    std::string sequence_dependence_model() const;
    void set_sequence_dependence_model(const std::string& model_name);
    void set_sequence_dependence_model(const StepParametersDB& steps_db,
                                       const ForceConstantsDB& fmat_db);

    // step sequence-dependent data accessors
    const StepSequence& bp_step_sequence(Size step_index) const;
    const TetramerSequence& bp_tetramer_sequence(Size step_index) const; //added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    const BpStepParams bp_step_intrinsic_parameters(Size step_index) const;
    const MatrixN& bp_step_force_constants(Size step_index) const;

    // update methods
    void update_collection_from_base_pairs(const std::vector<BasePair>& bps);
    void update_collection_from_bp_step_params(const std::vector<BpStepParams>&
                                               bp_step_params);
    void update_collection_from_bp_step_dofs(const std::vector<BpStepDofs>&
                                             bp_step_dofs);


private:

    // collection data
    std::vector<BasePair> m_base_pairs;
    std::vector<BpStepParams> m_bp_step_params;
    std::vector<BpStepDofs> m_bp_step_dofs;
    std::vector<SizePair> m_frozen_steps_domains;

    // sequence-dependence data
    std::vector<StepSequence> m_step_sequences;
    std::vector<TetramerSequence> m_tetramer_sequences; //added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    StepParametersDB m_step_seqdep;
    ForceConstantsDB m_fmat_seqdep;


};


#endif  // emDNA_BpCollection_h
