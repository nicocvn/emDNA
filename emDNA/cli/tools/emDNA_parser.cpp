// emDNA_parser tool
// Nicolas Clauvelin


#include <TCLAP_Includes.h>
#include <emDNA.h>
#include <emDNA_parser.h>


// entry point
int main(int argc, char* argv[]) {

    try {

        // command line parsing
        ParserData parser_data = parse_command_line(argc, argv);

        // bp collection
        BpCollection bp_collection;
        if (parser_data._input_x3DNAbp) {
            bp_collection =
            BpCollectionFactory::
            create_bp_collection_from_x3DNA_base_pairs_file(parser_data.
                                                            _filename);
//            bp_collection =
//            BpCollection_x3DNA::
//            create_bp_collection_from_base_pairs_file(parser_data._filename);
        }
        else if (parser_data._input_x3DNAparams) {
            bp_collection =
            BpCollectionFactory::
            create_bp_collection_from_x3DNA_bp_step_params_file(parser_data.
                                                                _filename);
//            bp_collection =
//            BpCollection_x3DNA::
//            create_bp_collection_from_step_parameters_file(parser_data.
//                                                           _filename);
        }
        else {

            // collection
            bp_collection =
            BpCollectionFactory::
            create_bp_collection_from_base_pairs_file(parser_data._filename);

            // sequence
            DS_ASSERT(parser_data._bp_list_sequence.size() ==
                      bp_collection.n_of_base_pairs(),
                      "incompatible sizes between BpCollection and sequence "
                      "string");
            bp_collection.set_collection_sequence(parser_data.
                                                  _bp_list_sequence);
        };

        // output
        if (parser_data._output_x3DNAbp) {
            std::stringstream output;
            BpCollection_x3DNA::format_as_base_pairs_file(output,
                                                          bp_collection);
            std::cout << output.str();
        }
        else if (parser_data._output_x3DNAparams) {
            std::stringstream output;
            BpCollection_x3DNA::format_as_step_parameters_file(output,
                                                               bp_collection);
            std::cout << output.str();
        }
        else {

            // loop over all base pairs
            for (Size i=0; i<bp_collection.n_of_base_pairs(); ++i) {
                std::cout << bp_collection.base_pair(i) << "\n";
            };
            
        };

        return 0;

    }

    catch(DNASim_ExitException& e) {
        exit(e._exit_code);
    };

};


// command line parsing option
ParserData parse_command_line(int argc, char* argv[]) {

    // command line parser
    TCLAP::CmdLine cmd_line("emDNA_parser - Nicolas Clauvelin, "
                            "Rutgers University",
                            '=',
                            "");

    // input arguments
    TCLAP::ValueArg<std::string> input_x3DNAbp("", "x3DNA-bp-input",
                                               "x3DNA base pairs input file.",
                                               true,
                                               "",
                                               "string");
    TCLAP::ValueArg<std::string> input_x3DNAparams("",
                                                   "x3DNA-bp-step-params-input",
                                                   "x3DNA bp step parameters "
                                                   "input file.",
                                                   true,
                                                   "",
                                                   "string");
    TCLAP::ValueArg<std::string> input_bplist("",
                                              "bp-list-input",
                                              "bplist "
                                              "input file.",
                                              true,
                                              "",
                                              "string");
    std::vector<TCLAP::Arg*> xor_input;
    xor_input.push_back(&input_x3DNAbp);
    xor_input.push_back(&input_x3DNAparams);
    xor_input.push_back(&input_bplist);
    cmd_line.xorAdd(xor_input);

    // output arguments
    TCLAP::SwitchArg output_x3DNAbp("", "get-x3DNA-bp",
                                    "x3DNA base pairs output.",
                                    true);
    TCLAP::SwitchArg output_x3DNAparams("", "get-x3DNA-params",
                                                "x3DNA bp step parameters "
                                                "output.",
                                                true);
    TCLAP::SwitchArg output_bplist("", "get-bp-list",
                                               "bp list output.",
                                               true);
    std::vector<TCLAP::Arg*> xor_output;
    xor_output.push_back(&output_x3DNAbp);
    xor_output.push_back(&output_x3DNAparams);
    xor_output.push_back(&output_bplist);
    cmd_line.xorAdd(xor_output);

    // sequence argument for bp list format
    TCLAP::ValueArg<std::string> bp_list_sequence("",
                                                  "bp-list-sequence",
                                                  "bplist sequence.",
                                                  false,
                                                  "",
                                                  "string");
    cmd_line.add(&bp_list_sequence);

    // parsing
    cmd_line.parse(argc, argv);

    // consistency check
    if (input_bplist.isSet() && !(bp_list_sequence.isSet())) {
        std::cout << "[[ missing sequence information ]]\n";
        std::cout << "using bplist format as input requires the specification ";
        std::cout << "of a sequence (--bp-list-sequence)\n";
        exit(-1);
    };

    // input settings
    ParserData parser_data;
    if (input_x3DNAbp.isSet()) {
        parser_data._input_x3DNAbp = true;
        parser_data._input_x3DNAparams = false;
        parser_data._input_bp_list = false;
        parser_data._filename = input_x3DNAbp.getValue();
    }
    else if (input_x3DNAparams.isSet()) {
        parser_data._input_x3DNAbp = false;
        parser_data._input_x3DNAparams = true;
        parser_data._input_bp_list = false;
        parser_data._filename = input_x3DNAparams.getValue();
    }
    else {
        parser_data._input_x3DNAbp = false;
        parser_data._input_x3DNAparams = false;
        parser_data._input_bp_list = true;
        parser_data._filename = input_bplist.getValue();
        parser_data._bp_list_sequence = bp_list_sequence.getValue();
    };

    // output settings
    if (output_x3DNAbp.isSet()) {
        parser_data._output_x3DNAbp = true;
        parser_data._output_x3DNAparams = false;
        parser_data._output_bp_list = false;
    }
    else if (output_x3DNAparams.isSet()) {
        parser_data._output_x3DNAbp = false;
        parser_data._output_x3DNAparams = true;
        parser_data._output_bp_list = false;
    }
    else {
        parser_data._output_x3DNAbp = false;
        parser_data._output_x3DNAparams = false;
        parser_data._output_bp_list = true;
    };

    return parser_data;

};
