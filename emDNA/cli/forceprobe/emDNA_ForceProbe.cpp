// emDNA_ForceProbe app
// Nicolas Clauvelin


#include <emDNA.h>
#include <ForceRampMinimizer.h>
#include <emDNA_ForceProbe.h>


// entry point
int main(int argc, char* argv[]) {

    try {

        // input parsing
        AppData input_data = parse_input(argc, argv);

        // base collection
        BpCollection base_coll = create_bp_collection(input_data._base_input);
        base_coll.set_sequence_dependence_model(input_data._base_seqdep_model);
        base_coll.
        set_frozen_steps_domains(emDNA_Utils::
                                 parse_frozen_steps_vector(input_data.
                                                           _base_bound_domains)
                                 );

        // force ramp
        ForceRamp force_ramp;
        Vector3 force_probe_data(input_data._force_probe);
        force_ramp._initial_force = force_probe_data[X];
        force_ramp._final_force = force_probe_data[Y];
        force_ramp._n_increments = force_probe_data[Z];

        // normalized force direction
        force_ramp._force_direction = Vector3(input_data._force_probe_dir);
        force_ramp._force_direction.normalize();

        // message
        print_inspector_output(base_coll, force_ramp);
        std::cout << "\n";

        // minimization
        ForceRampMinimizer minimizer;
        minimizer.set_initial_bp_collection(base_coll);
        minimizer.set_force_ramp(force_ramp);
        minimizer.force_ramp_minimization();

        return 0;

    }
    catch (DNASim_ExitException& e) {
        exit(e._exit_code);
    };

};


// input parsing function
AppData parse_input(int argc, char* argv[]) {

    if (argc != 2) {
        std::cout << "wrong syntax\n";
        std::cout << "usage: emDNA_ForceProbe input_file\n";
        throw DNASim_ExitException(-1);
    };

    // input file reading
    std::map<std::string,std::string> input_data =
    InputFileHandler::read_key_value_file(argv[1], ".", '=', '#');

    // input data
    AppData app_data;

    // base_collection
    if (!find_key_value(input_data,
                        "base_collection",
                        app_data._base_input)) {
        std::cout << "base_configuration key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    // base_bound_domains
    if (!find_key_value(input_data,
                        "base_bound_domains",
                        app_data._base_bound_domains)) {
        app_data._base_bound_domains = std::string();
    };

    // base_seqdep_model
    if (!find_key_value(input_data,
                        "base_seqdep_model",
                        app_data._base_seqdep_model)) {
        app_data._base_seqdep_model = "IdealDNA";
    };

    // force_probe
    if (!find_key_value(input_data,
                        "force_probe",
                        app_data._force_probe)) {
        std::cout << "force_probe key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    // force_probe_direction
    if (!find_key_value(input_data,
                        "force_probe_direction",
                        app_data._force_probe_dir)) {
        std::cout << "force_probe_direction key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    return app_data;

};


// bp collection creation function
BpCollection create_bp_collection(const std::string& base_input) {

    // tokenizing
    // first token is input  mode and second one is filename
    std::vector<std::string> tokens =
    EnhancedString::tokenize_string(base_input, ':');
    DS_ASSERT(tokens.size()==2,
              "wrong format for base_collection key\n");

    // x3DNAbp input mode
    if (tokens[0] == "x3DNAbp") {
        return
        BpCollectionFactory::
        create_bp_collection_from_x3DNA_base_pairs_file(tokens[1]);
    }

    // x3DNAparams input mode
    else if (tokens[0] == "x3DNAparams") {
        return
        BpCollectionFactory::
        create_bp_collection_from_x3DNA_bp_step_params_file(tokens[1]);
    }

    // bp list input mode
    else if (tokens[0] == "bplist") {
        BpCollection bp_collection =
        BpCollectionFactory::
        create_bp_collection_from_base_pairs_file(tokens[1]);
        bp_collection.set_collection_dummy_sequence();
        return bp_collection;
    }
    else {
        std::cout << "unsupported input mode for base_collection\n";
        throw DNASim_ExitException(-1);
    };
    
};


// inspector message function
void print_inspector_output(const BpCollection& base_collection,
                            const ForceRamp& force_ramp) {

    // base collection size
    const Size& base_n_steps = base_collection.n_of_bp_steps();

    // base sequence dependence model
    const std::string& base_seqdep =
    base_collection.sequence_dependence_model();

    // base frozen domains
    const std::vector<SizePair>& frozen_domains =
    base_collection.frozen_steps_domains();

    // base collection message
    std::cout << "--- base bp collection:\n";
    std::cout << "  number of steps: " << base_n_steps << "\n";
    std::cout << "  sequence-dependence model: " << base_seqdep << "\n";
    std::cout << "  bound domains: ";
    print_vector(frozen_domains);
    std::cout << "\n";

    // binding message
    std::cout << "--- force ramp:\n";
    std::cout << "  pulling direction: " << force_ramp._force_direction << "\n";
    std::cout << "  initial force: " << force_ramp._initial_force << "\n";
    std::cout << "  final force: " << force_ramp._final_force << "\n";
    std::cout << "  # increments: " << force_ramp._n_increments << "\n";
    std::cout << "\n";
    
};

