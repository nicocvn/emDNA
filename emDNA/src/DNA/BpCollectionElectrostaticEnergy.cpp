// BpCollectionElectrostaticEnergy class
// Juan Wei, Nicolas Clauvelin


#include "DNA/DNAElectrostaticsParams.h"
#include "DNA/BpCollection.h"
#include "DNA/BpCollectionElectrostaticEnergy.h"


Real
BpCollectionElectrostaticEnergy::electrostatic_energy(const BpCollection&
                                                      bp_collection) {
    
    Real E(FLOAT_INIT);

    // list of base pairs and iterators
	const std::vector<BasePair>& bps = bp_collection.base_pairs();
    const Size n_bps = bps.size();
	
	// loop over all base pairs
	// July 2019: include a "cheat" for circular DNA (lines 22-25; change i=0 to i=1 in line 26
	for (Size j=DNAElec::ExclusionRange; j<n_bps-1; ++j) {
       const Real rij = (bps[j].origin()-bps[0].origin()).norm();
       E += std::exp(-DNAElec::kappaDebye*rij)/rij;
    };
    for (Size i=1; i<n_bps; ++i) {
        for (Size j=i+DNAElec::ExclusionRange; j<n_bps; ++j) {
            const Real rij = (bps[j].origin()-bps[i].origin()).norm();
            E += std::exp(-DNAElec::kappaDebye*rij)/rij;
        };
    };
		
	E *= DNAElec::dhConstant;

    return E;

};
