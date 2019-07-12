// BpCollectionFactory class
// Nicolas Clauvelin


#include <BpCollection.h>
#include <BpCollection_x3DNA.h>
#include <BpCollectionFactory.h>


// static factory method - base pairs file
BpCollection BpCollectionFactory::
create_bp_collection_from_base_pairs_file(const std::string& filename) {

    std::vector<BasePair> base_pairs(read_input_file<BasePair>(filename));

    // orthogonalization
    const Size n_bp = base_pairs.size();
    for (Size i=0; i<n_bp; i++)
        base_pairs[i].orthogonalize();

    // new collection with dummy sequence
    // bplist format does not carry sequence information
    BpCollection bp_coll(BpCollection::collection_from_base_pairs(base_pairs));
    bp_coll.set_collection_dummy_sequence();

    return bp_coll;

};


// static factory method - x3DNA base pairs file
BpCollection BpCollectionFactory::
create_bp_collection_from_x3DNA_base_pairs_file(const std::string& filename) {
    return
    BpCollection_x3DNA::create_bp_collection_from_base_pairs_file(filename);
};


// static factory method - x3DNA bp step params file
BpCollection BpCollectionFactory::
create_bp_collection_from_x3DNA_bp_step_params_file(const std::string&
                                                    filename) {
    return BpCollection_x3DNA::
    create_bp_collection_from_step_parameters_file(filename);
};


// file reading static method
template <class DataType>
std::vector<DataType> BpCollectionFactory::read_input_file(const std::string&
                                                           filename) {

    // file reading
    std::vector<std::string> file_lines(InputFileHandler::read_file(filename,
                                                                    "."));

    // data casting
    std::vector<DataType> data_vector;
    data_vector.reserve(file_lines.size());
    for (Size i=0; i<file_lines.size(); ++i)
        if (file_lines[i].size() != 0 && file_lines[i][0] != '#') {
            data_vector.push_back(DataType(file_lines[i]));
        };

    return data_vector;

};
