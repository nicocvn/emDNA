// BpCollectionElasticEnergy class
// Nicolas Clauvelin


#include <BpCollection.h>
#include <IndexManager.h>
#include <BpCollectionElasticEnergy.h>


// bp collection elastic energy function
// the elastic energy is computed over the free bp steps
Real
BpCollectionElasticEnergy::elastic_energy(const BpCollection& bp_collection,
                                          const IndexManager& index_manager) {

    // energy
    Real E(FLOAT_INIT);

    // loop over free steps
    auto begin = index_manager.free_step_indexes_begin();
    auto end = index_manager.free_step_indexes_end();
    for (auto free_it = begin; free_it != end; ++free_it) {

        E +=
        BpCollectionElasticEnergy::
        single_step_elastic_energy(free_it->collection_index(),
                                   bp_collection);

    };

    return E;

};


// single step elastic energy
Real
BpCollectionElasticEnergy::single_step_elastic_energy(Size step_index,
                                                      const BpCollection&
                                                      bp_collection) {

    const VectorN dp(bp_collection.bp_step_params(step_index).inline_vector()
                     -bp_collection.bp_step_intrinsic_parameters(step_index).
                     inline_vector());

    return Real(0.5)*dp.dot(bp_collection.
                            bp_step_force_constants(step_index)*dp);

};


// bp collection elastic energy splits
MatrixN
BpCollectionElasticEnergy::elastic_energy_splits(const BpCollection&
                                                 bp_collection,
                                                 const IndexManager&
                                                 index_manager) {

    // container
    MatrixN energy_split(emDNAConstants::StepParametersDim);

    // loop over the free steps to compute the splits
    const std::vector<Size> free_steps =
    index_manager.free_step_collection_indices();
    for (const Size& idx : free_steps) {

        // intrinsic step parameters
        const BpStepParams& p0 =
        bp_collection.bp_step_intrinsic_parameters(idx);

        // force constants
        const MatrixN& fmat =
        bp_collection.bp_step_force_constants(idx);

        // inline vectors
        const VectorN pvec =
        bp_collection.bp_step_params(idx).inline_vector();
        const VectorN pvec0 = p0.inline_vector();

        // loop over contributions
        for (Size i=0; i<emDNAConstants::StepParametersDim; ++i) {
            for (Size j=0; j<emDNAConstants::StepParametersDim; ++j) {

                energy_split(i,j) +=
                Real(0.5)*fmat(i,j)*(pvec[i]-pvec0[j])*(pvec[i]-pvec0[j]);

            };
        };
        
    };
    
    return energy_split;


};

