// IndexManager class
// Nicolas Clauvelin


#include "DNA/IndexManager.h"


// class defautl constructor
IndexManager::IndexManager() :
m_n_steps(),
m_free_steps(),
m_frozen_steps() {};


// class constructor by copy
IndexManager::IndexManager(const IndexManager& idx_mgr) :
m_n_steps(idx_mgr.m_n_steps),
m_free_steps(idx_mgr.m_free_steps),
m_frozen_steps(idx_mgr.m_frozen_steps) {};


// class destructor
IndexManager::~IndexManager() {};


// copy operator
IndexManager& IndexManager::operator=(const IndexManager& idx_mgr) {
    m_n_steps = idx_mgr.m_n_steps;
    m_free_steps = idx_mgr.m_free_steps;
    m_frozen_steps = idx_mgr.m_frozen_steps;
    return *this;
};


// total size setup method
void IndexManager::set_n_of_bp_step(Size n_bp_steps) {
    m_n_steps = n_bp_steps;
};


// frozen step indexes setup method
void IndexManager::set_frozen_steps_indexes(const std::vector<SizePair>&
                                            frozen_step_domains) {

    // cleaning
    m_free_steps.clear();
    m_frozen_steps.clear();

    // counters
    Size n_free = 0;
    Size n_frozen = 0;

    // loop over all steps
    for (Size i=0; i<m_n_steps; ++i) {

        // check if the index is in a frozen domain
        bool is_frozen = false;
        for (const SizePair& domain : frozen_step_domains) {

            // this step is frozen
            // store the index as frozen, increment the counter and break
            if ((domain.first <= i) && (i <= domain.second)) {
                m_frozen_steps.push_back(IndexPair(i, n_frozen));
                ++n_frozen;
                is_frozen = true;
                break;
            };

        };

        // non-frozen step
        // store the index as free and increment the counter
        if (is_frozen == false) {
            m_free_steps.push_back(IndexPair(i, n_free));
            ++n_free;
        };

    };

    // intregrity check
    DS_ASSERT(m_free_steps.size() + m_frozen_steps.size() == m_n_steps,
              "index error in step index manager setup method");
    DS_ASSERT(m_free_steps.size() == n_free,
              "index error in step index manager setup method");
    DS_ASSERT(m_frozen_steps.size() == n_frozen,
              "index error in step index manager setup method");

};


// size accessors
Size IndexManager::n_of_dofs_variables() const {
    return m_n_steps*emDNAConstants::StepParametersDim;
};
Size IndexManager::n_of_free_dofs_variables() const {
    return m_free_steps.size()*emDNAConstants::StepParametersDim;
};
Size IndexManager::n_of_parameters_variables() const {
    return m_n_steps*emDNAConstants::StepParametersDim;
};
Size IndexManager::n_of_free_parameters_variables() const {
    return m_free_steps.size()*emDNAConstants::StepParametersDim;
};


// number of free step accessor
Size IndexManager::n_of_free_steps() const {
    return m_free_steps.size();
};


// all free dofs indices accessor
std::vector<Size> IndexManager::inline_free_dofs_indices() const {

    // number of free dofs variables
    std::vector<Size> inline_indexes(n_of_free_dofs_variables());

    // iterators
    const auto free_begin = free_step_indexes_begin();
    const auto free_end = free_step_indexes_end();

    // indices
    for (auto free_it = free_begin; free_it != free_end; ++free_it) {
        const Size coll_index = free_it->collection_index();
        const Size local_index = free_it->local_index();
        for (Size i=0; i<emDNAConstants::StepParametersDim; ++i)
            inline_indexes[local_index*emDNAConstants::StepParametersDim + i] =
            coll_index*emDNAConstants::StepParametersDim + i;
    };

    return inline_indexes;

};


// all free steps collection indices accessor
std::vector<Size> IndexManager::free_step_collection_indices() const {
    std::vector<Size> indices;
    indices.reserve(m_free_steps.size());
    std::for_each(m_free_steps.begin(), m_free_steps.end(),
                  [&indices](const IndexPair& p) {
                      indices.push_back(p.collection_index());
                  });
    return indices;
};


// free steps iterators
IndexPairConstIt IndexManager::free_step_indexes_begin() const {
    return m_free_steps.begin();
};
IndexPairConstIt IndexManager::free_step_indexes_end() const {
    return m_free_steps.end();
};


// bc reduced step index accessor
// the bc reduced step is simply an alias to the last free step
const IndexPair& IndexManager::bc_reduced_step() const {
    return m_free_steps.back();
};


// number of frozen step accessors
Size IndexManager::n_of_frozen_steps() const {
    return m_frozen_steps.size();
};


// all frozen steps collection indices accessor
std::vector<Size> IndexManager::frozen_step_collection_indices() const {
    std::vector<Size> indices;
    indices.reserve(m_frozen_steps.size());
    std::for_each(m_frozen_steps.begin(), m_frozen_steps.end(),
                  [&indices](const IndexPair& p) {
                      indices.push_back(p.collection_index());
                  });
    return indices;
};


// free steps iterators
IndexPairConstIt IndexManager::frozen_step_indexes_begin() const {
    return m_frozen_steps.begin();
};
IndexPairConstIt IndexManager::frozen_step_indexes_end() const {
    return m_frozen_steps.end();
};


// dof and parameter global coordinate methods
Size IndexManager::dof_global_coordinate(Size step_index,
                                         BpStepDof dof_index) const {
    return (Size)dof_index+step_index*emDNAConstants::StepParametersDim;
};
Size IndexManager::parameter_global_coordinate(Size step_index,
                                               BpStepParameter param_index)
const {
    return (Size)param_index+step_index*emDNAConstants::StepParametersDim;
};

