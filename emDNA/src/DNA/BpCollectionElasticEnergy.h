// BpCollectionElasticEnergy class
// Nicolas Clauvelin


#ifndef emDNA_BpCollectionElasticEnergy_h
#define emDNA_BpCollectionElasticEnergy_h


#include <emDNA_Includes.h>
class BpCollection;
class IndexManager;


class BpCollectionElasticEnergy {


public:

    // bp collection elastic energy function
    static Real elastic_energy(const BpCollection& bp_collection,
                               const IndexManager& index_manager);

    // single step elastic energy
    static Real single_step_elastic_energy(Size step_index,
                                           const BpCollection& bp_collection);

    // bp collection elastic energy splits
    static MatrixN elastic_energy_splits(const BpCollection& bp_collection,
                                         const IndexManager& index_manager);


private:

    // private constructors and copy operator
    // this class is not designed to be instantiated
    BpCollectionElasticEnergy() = delete;
    BpCollectionElasticEnergy(const BpCollectionElasticEnergy&
                              bp_coll_eng) = delete;
    BpCollectionElasticEnergy(BpCollectionElasticEnergy&&
                              bp_coll_eng) = delete;
    BpCollectionElasticEnergy& operator=(const BpCollectionElasticEnergy&
                                         bp_coll_eng) = delete;
    BpCollectionElasticEnergy& operator=(BpCollectionElasticEnergy&&
                                         bp_coll_eng) = delete;
    ~BpCollectionElasticEnergy() = delete;


};


#endif  // emDNA_BpCollectionElasticEnergy_h
