// emProBind app
// Nicolas Clauvelin


#include <emDNA_Utils.h>
#include <RampMinimizerEEDR.h>
#include <RampMinimizerPulling.h>
#include <BpCollection.h>
#include <BpCollectionFactory.h>
#include <emProBind.h>


// entry point
int main(int argc, char* argv[]) {

    try {

        // input parsing
        emProBindData input_data = parse_input(argc, argv);

        // base collection
        BpCollection base_coll = create_bp_collection(input_data._base_input);
        base_coll.set_sequence_dependence_model(input_data._base_seqdep_model);
        base_coll.
        set_frozen_steps_domains(emDNA_Utils::
                                 parse_frozen_steps_vector(input_data.
                                                           _base_bound_domains)
                                 );

        // protein collection
        BpCollection protein_coll = create_bp_collection(input_data.
                                                         _protein_input);

        // binding sites
        std::vector<Size> sites =
        parse_protein_binding_sites(input_data._protein_binding_sites);

        // message
        print_inspector_output(base_coll, protein_coll, sites);
        std::cout << "[[ ramp binding with " << input_data._ramp_sampling;
        std::cout << " steps ]]\n\n";

        // collection type message
        if (input_data._coll_type == CollectionType::ClampedCollection)
            std::cout << "[[ EEDR minimization ]]\n\n";
        if (input_data._coll_type == CollectionType::PullCollection) {
            std::cout << "[[ pulling force minimization ]]\n";
            std::cout << "[[ pulling force = " << input_data._pulling_force
                << " ]]\n\n";
        };
        std::cout << "\n";

        // binding ramp minimizer instanciation
        std::unique_ptr<RampMinimizer> ramp_minim_ptr = nullptr;
        if (input_data._coll_type == CollectionType::ClampedCollection) {
            ramp_minim_ptr =
            std::unique_ptr<RampMinimizer>(new RampMinimizerEEDR());
        };
        if (input_data._coll_type == CollectionType::PullCollection) {
            Vector3 f(input_data._pulling_force);
            ramp_minim_ptr =
            std::unique_ptr<RampMinimizerPulling>(new RampMinimizerPulling(f));
        };

        // minimization
        ramp_minim_ptr->set_initial_bp_collection(base_coll);
        ramp_minim_ptr->set_protein_bp_collection(protein_coll);
        ramp_minim_ptr->set_binding_ramp_sampling(input_data._ramp_sampling);
        ramp_minim_ptr->set_protein_binding_indices(sites);
        bool flag = ramp_minim_ptr->binding_ramp_minimization();

        // message
        if (flag)
            std::cout << "binding ramp minimization successful\n";
        else
            std::cout << "!!! binding ramp minimization not successful !!!\n";

        return 0;

    }

    catch (DNASim_ExitException& e) {
        exit(e._exit_code);
    };
    
};


// input parsing function
emProBindData parse_input(int argc, char* argv[]) {

    if (argc != 2) {
        std::cout << "wrong syntax\n";
        std::cout << "usage: emDNA_probind input_file\n";
        throw DNASim_ExitException(-1);
    };

    // input file reading
    std::map<std::string,std::string> input_data =
    InputFileHandler::read_key_value_file(argv[1], ".", '=', '#');

    // input data
    emProBindData pro_binding_data;

    // collection type
    std::string coll_type;
    if (!find_key_value(input_data,
                        "collection_type",
                        coll_type)) {
        std::cout << "collection_tupe key missing in input file\n";
        throw DNASim_ExitException(-1);
    };
    if (coll_type == "EEDR")
        pro_binding_data._coll_type = CollectionType::ClampedCollection;
    else if (coll_type == "pulling")
        pro_binding_data._coll_type = CollectionType::PullCollection;
    else
        DS_ASSERT(false, "unknown collection type in input file");

    // pulling force
    if (!find_key_value(input_data,
                        "pulling_force",
                        pro_binding_data._pulling_force)) {
        pro_binding_data._pulling_force = std::string();
    };

    // base_collection
    if (!find_key_value(input_data,
                        "base_collection",
                        pro_binding_data._base_input)) {
        std::cout << "base_configuration key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    // base_bound_domains
    if (!find_key_value(input_data,
                        "base_bound_domains",
                        pro_binding_data._base_bound_domains)) {
        pro_binding_data._base_bound_domains = std::string();
    };

    // base_seqdep_model
    if (!find_key_value(input_data,
                        "base_seqdep_model",
                        pro_binding_data._base_seqdep_model)) {
        pro_binding_data._base_seqdep_model = "IdealDNA";
    };

    // protein_collection
    if (!find_key_value(input_data,
                        "protein_collection",
                        pro_binding_data._protein_input)) {
        std::cout << "protein_collection key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    // protein_binding_sites
    if (!find_key_value(input_data,
                        "protein_binding_sites",
                        pro_binding_data._protein_binding_sites)) {
        std::cout << "protein_binding_sites key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    // binding_ramp_sampling
    if (!find_key_value(input_data,
                        "binding_ramp_sampling",
                        pro_binding_data._ramp_sampling)) {
        std::cout << "binding_ramp_sampling key missing in input file\n";
        throw DNASim_ExitException(-1);
    };

    return pro_binding_data;
    
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
                            const BpCollection& protein_collection,
                            const std::vector<Size>& binding_indices) {

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

    // protein collection size
    const Size& protein_n_steps = protein_collection.n_of_bp_steps();

    // protein collection message
    std::cout << "--- protein bp collection:\n";
    std::cout << "  number of steps: " << protein_n_steps << "\n";

    // binding message
    std::cout << "--- binding information:\n";
    std::cout << "  # proteins: " << binding_indices.size() << "\n";
    std::cout << "  binding sites: ";
    print_vector(binding_indices);
    std::cout << "\n\n";

};


// protein_binding_sites key parser function
std::vector<Size> parse_protein_binding_sites(const std::string& s) {

    // the string is parse a vector
    VectorN real_sites(s);

    // casting
    std::vector<Size> sites;
    sites.reserve(real_sites.size());
    for (Size i=0, end=real_sites.size(); i<end; ++i)
        sites.push_back(static_cast<Size>(real_sites[i]-Real(1)));

    return sites;
    
};
