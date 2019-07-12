// BpCollection_x3DNA class
// Nicolas Clauvelin


#include <BpCollection.h>
#include <BpCollection_x3DNA.h>


// x3DNA reading method for step parameters file
BpCollection BpCollection_x3DNA::
create_bp_collection_from_step_parameters_file(const std::string& filename) {

    // parsed data
    Data_x3DNA<StepParameters> parsed_data =
    Parser_x3DNA::read_step_parameters_file(filename);

    // convert StepParameters to BpStepParams
    const std::vector<StepParameters>& prms_vec = parsed_data._data;
    std::vector<BpStepParams> convert_prms;
    convert_prms.reserve(prms_vec.size());
    for (auto it = prms_vec.begin(); it != prms_vec.end(); ++it)
        convert_prms.push_back(BpStepParams(it->inline_vector()));

    // bp collection
    BpCollection bp_collection =
    BpCollection::collection_from_bp_step_params(convert_prms, BasePair());
    bp_collection.set_collection_sequence(parsed_data._sequence);

    return bp_collection;

};


// x3DNA reading method for base pairs file
BpCollection BpCollection_x3DNA::
create_bp_collection_from_base_pairs_file(const std::string& filename) {

    // parsed data
    Data_x3DNA<Triad> parsed_data =
    Parser_x3DNA::read_base_pairs_file(filename);

    // bp collection
    BpCollection bp_collection =
    BpCollection::collection_from_base_pairs(parsed_data._data);
    bp_collection.set_collection_sequence(parsed_data._sequence);

    return bp_collection;

};


// x3DNA output method for step parameters format
void BpCollection_x3DNA::
format_as_step_parameters_file(std::stringstream& output,
                               const BpCollection& bp_collection) {

    // number of steps
    const Size n = bp_collection.n_of_bp_steps();

    // convert BpStepParams to StepParameters
    std::vector<StepParameters> prms;
    prms.reserve(n);
    for (Size i=0; i<n; ++i)
        prms.
        push_back(StepParameters(bp_collection.
                                 bp_step_params(i).inline_vector()));

    // output data
    Data_x3DNA<StepParameters> output_data;
    output_data._data = prms;
    output_data._sequence = bp_collection.collection_sequence();

    // output
    Parser_x3DNA::create_step_parameters_output(output, output_data);

};


// x3DNA output method for base pairs format
void BpCollection_x3DNA::format_as_base_pairs_file(std::stringstream& output,
                                                   const BpCollection&
                                                   bp_collection) {

    // output data
    Data_x3DNA<Triad> output_data;
    output_data._data = bp_collection.base_pairs();
    output_data._sequence = bp_collection.collection_sequence();

    // output
    Parser_x3DNA::create_base_pairs_output(output, output_data);

};


// x3DNA writing methods
void BpCollection_x3DNA::
write_as_step_parameters_file(const std::string& filename,
                              const BpCollection& bp_collection) {

    // output file
    OutputFileHandler output_file(filename, ".");

    // formatting
    std::stringstream output;
    format_as_step_parameters_file(output, bp_collection);

    // file writing
    output_file.open();
    output_file.write_stream(output);
    output_file.close();

};
void BpCollection_x3DNA::
write_as_base_pairs_file(const std::string& filename,
                         const BpCollection& bp_collection) {

    // output file
    OutputFileHandler output_file(filename, ".");

    // formatting
    std::stringstream output;
    format_as_base_pairs_file(output, bp_collection);

    // file writing
    output_file.open();
    output_file.write_stream(output);
    output_file.close();

};

