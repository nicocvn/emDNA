// BpCollectionElectrostaticEnergy class
// Juan Wei, Nicolas Clauvelin


#ifndef emDNA_BpCollectionElectrostaticEnergy_h
#define emDNA_BpCollectionElectrostaticEnergy_h


#include <emDNA_Includes.h>


class BpCollection;


class BpCollectionElectrostaticEnergy {


public:

	static Real electrostatic_energy(const BpCollection& bp_collection);


};


#endif  // emDNA_BpCollectionElectrostaticEnergy_h
