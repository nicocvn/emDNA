// StepGradient structure
// Nicolas Clauvelin


#include <IndexManager.h>
#include <StepGradient.h>


// flattening function for StepGradientVec
VectorN flatten_StepGradientVec(const StepGradientVec& v) {

    // stl container
    std::vector<Real> vec;
    vec.reserve(v.size()*emDNAConstants::StepParametersDim);

    // flattening
    std::for_each(v.begin(), v.end(), [&vec](const StepGradient& g){
        vec.insert(vec.end(),
                   g._rotation.packed_values()._values,
                   g._rotation.packed_values()._values+3);
        vec.insert(vec.end(),
                   g._translation.packed_values()._values,
                   g._translation.packed_values()._values+3);
    });

    // VectorN container
    return VectorN(vec);

};


// filtering function for a vector of StepGradient
StepGradientVec filter_free_dofs_grads(const StepGradientVec& v,
                                       const IndexManager& idx_mgr) {

    StepGradientVec filtered_grads;
    filtered_grads.reserve(idx_mgr.n_of_free_steps());

    const auto free_begin = idx_mgr.free_step_indexes_begin();
    const auto free_end = idx_mgr.free_step_indexes_end();
    for (auto it = free_begin; it != free_end; ++it) {
        filtered_grads.push_back(v[it->collection_index()]);
    };

    return filtered_grads;

};
