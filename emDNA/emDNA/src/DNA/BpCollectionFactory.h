// BpCollectionFactory class
// Nicolas Clauvelin


#ifndef emDNA_BpCollectionFactory_h
#define emDNA_BpCollectionFactory_h


#include <emDNA_Includes.h>
class BpCollection;


class BpCollectionFactory {


public:

    // static factory methods
    static BpCollection
    create_bp_collection_from_base_pairs_file(const std::string& filename);
    static BpCollection
    create_bp_collection_from_x3DNA_base_pairs_file(const std::string&
                                                    filename);
    static BpCollection
    create_bp_collection_from_x3DNA_bp_step_params_file(const std::string&
                                                        filename);


private:

    // file reading method
    template <class DataType>
    static std::vector<DataType> read_input_file(const std::string& filename);

    // private constructors and copy operator
    BpCollectionFactory() = delete;
    BpCollectionFactory(const BpCollectionFactory&
                        bp_collection_factory) = delete;
    BpCollectionFactory(BpCollectionFactory&&
                        bp_collection_factory) = delete;
    BpCollectionFactory& operator=(const BpCollectionFactory&
                                   bp_collection_factory) = delete;
    BpCollectionFactory& operator=(BpCollectionFactory&&
                                   bp_collection_factory) = delete;
    ~BpCollectionFactory() = delete;
    


};


#endif  // emDNA_BpCollectionFactory_h
