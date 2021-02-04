// BpCollectionElasticEnergy unit tests
// Nicolas Clauvelin


#include <emDNA.h>
#include "utest_BpCollectionElasticEnergy.h"


namespace {


    // ElasticEnergyComputation
    TEST_F(BpCollectionElasticEnergyTest, ElasticEnergyComputation) {

        BpCollection bpcoll =
        BpCollection::collection_from_base_pairs(test_base_pairs);
        bpcoll.set_collection_dummy_sequence();
        bpcoll.set_sequence_dependence_model("AnisoDNA_304");

        IndexManager idx_mgr;
        idx_mgr.set_n_of_bp_step(bpcoll.n_of_bp_steps());
        idx_mgr.set_frozen_steps_indexes(std::vector<SizePair>());

        Real E = BpCollectionElasticEnergy::elastic_energy(bpcoll, idx_mgr);

        EXPECT_NEAR(exactE, E, ZERO_TOL);

    };

}
