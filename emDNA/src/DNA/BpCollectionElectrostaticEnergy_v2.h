// BpCollectionElectrostaticEnergy class
// Juan Wei, Nicolas Clauvelin


#include <DNAElectrostaticsParams.h>
#include <BpCollection.h>
#include <BpCollectionElectrostaticEnergy.h>

Real
BpCollectionElectrostaticEnergy::electrostatic_energy(const BpCollection&
                                                      bp_collection) {

    Real E(FLOAT_INIT);

    // list of base pairs and iterators
    const std::vector<BasePair>& bps = bp_collection.base_pairs();
    const Size n_bps = bps.size();

    // loop over all base pairs
    //for (Size j=DNAElec::ExclusionRange; j<n_bps; ++j)
    //{
    //    const Real rij = (bps[j].origin()-bps[0].origin()).norm();
    //    E += std::exp(-DNAElec::kappaDebye*rij)/rij;
    //};


    // ---NEW--- 11 Nov 2020
    // check if current bpcollection is circular
    bool circular = false;

    if (  ) {
        circular = true;
    };

    // iterate over all base pairs
    if (circular == false) {
        for (Size i=0; i<n_bps; ++i) {
            for (Size j=i+DNAElec::ExclusionRange; j<n_bps; ++j) {
                const Real rij = (bps[j].origin()-bps[i].origin()).norm();
                E += std::exp(-DNAElec::kappaDebye*rij)/rij;
            };
        };
    } else if (circular == true) {
        for (Size i=0; i<n_bps; ++i) {
            const int x = std::abs(i-DNAElec::ExclusionRange);

            for (Size j=i+DNAElec::ExclusionRange; j<n_bps; ++j) {
                if (j <= n_bps-x) {
                    const Real rij = (bps[j].origin()-bps[i].origin()).norm();
                    E += std::exp(-DNAElec::kappaDebye*rij)/rij;
                };
            };
        };
    }

    E *= DNAElec::dhConstant;
    return E;

};
