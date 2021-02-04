// emDNA_Application class
// Nicolas Clauvelin


#include <emDNA_CommandLine.h>


#include <emDNA.h>
#include <emDNA_CommandLine.h>
#include <emDNA_Application.h>

// bp collection interfaces
//#include <FreeBpCollection.h>
//#include <AnchoredBpCollection.h>
//#include <ClampedBpCollection.h>
//#include <PullingBpCollection.h>

// needed for external force field

// needed for exception
//#include <Alglib_Includes.h>


// class default constructor
emDNA_Application::emDNA_Application() {};


// class destructor
emDNA_Application::~emDNA_Application() {};


// application command line parsing
void emDNA_Application::parse_command_line(int argc, char* argv[]) {
    m_cmd_line.parse_command_line_input(argc, argv);
};


// application setup method
void emDNA_Application::setup_application() {

    // log file
    setup_log_file(m_cmd_line.output_name());

    // log command line input
    m_log_file.log_emDNA_CommandLine(m_cmd_line);

    // bp collection
    setup_bp_collection();

    // log bp collection
    m_log_file.log_BpCollection(m_bp_collection);

    // setup message
    if (!m_cmd_line.quiet_flag())
        application_setup_message();

};


// minimization method
void emDNA_Application::perform_minimization() {

    // callback options
    CallbackOptions callback_opts;
    callback_opts._callback_output = !(m_cmd_line.quiet_flag());
    callback_opts._output_current_energy = m_cmd_line.energy_progress_flag();
    callback_opts._output_collection_frequency =
    m_cmd_line.output_progress_freq();

    // minim settings
    AlglibMinSettings minim_settings;
    AlgligMinSettings_Ptr minim_settings_ptr = nullptr;
    if (m_cmd_line.minim_settings(minim_settings))
        minim_settings_ptr =
        AlgligMinSettings_Ptr(new AlglibMinSettings(minim_settings));

    // message
    if (!m_cmd_line.quiet_flag())
        minimization_pre_flight_message(minim_settings_ptr);

    // minimization
    MinimizationResults minim_res =
    MinimizerAgent::alglib_minimization(m_bp_intf_ptr,
                                        callback_opts,
                                        minim_settings_ptr);

    // post-processing
    m_bp_collection = minim_res._optimized_bp_collection;

    // log minim results
    m_log_file.log_MinimResults(minim_res);
    m_log_file.log_EnergyComposition(minim_res);

    // minimization message
    if (!m_cmd_line.quiet_flag())
        minimization_post_flight_message(minim_res);

};


// results output method
void emDNA_Application::output_optimized_bp_collection() {

    // output file name
    std::string opt_filename = m_cmd_line.output_name() + "_opt.txt";

    // formatting
    if (m_cmd_line.input_mode() == InputMode::x3DNA_Bp) {
        BpCollection_x3DNA::write_as_base_pairs_file(opt_filename,
                                                     m_bp_collection);
    };
    if (m_cmd_line.input_mode() == InputMode::x3DNA_Params) {
        BpCollection_x3DNA::write_as_step_parameters_file(opt_filename,
                                                          m_bp_collection);
    };
    if (m_cmd_line.input_mode() == InputMode::BpList) {

        // output file
        OutputFileHandler output(opt_filename, ".");

        // loop over all base pairs
        output.open();
        for (Size i=0; i<m_bp_collection.n_of_base_pairs(); ++i) {
            output.stream() << m_bp_collection.base_pair(i) << "\n";
        };
        output.close();

    };

};


// bp collection setup method
void emDNA_Application::setup_bp_collection() {

    // input filename
    std::string input_filename = m_cmd_line.input_file();

    // bp collection from x3DNA bp input file
    if (m_cmd_line.input_mode() == InputMode::x3DNA_Bp)
        m_bp_collection =
        BpCollectionFactory::
        create_bp_collection_from_x3DNA_base_pairs_file(input_filename);

    // bp collection from x3DNA step params file
    else if (m_cmd_line.input_mode() == InputMode::x3DNA_Params)
        m_bp_collection =
        BpCollectionFactory::
        create_bp_collection_from_x3DNA_bp_step_params_file(input_filename);

    // bp collection from bp list input file
    else if (m_cmd_line.input_mode() == InputMode::BpList)
        m_bp_collection =
        BpCollectionFactory::
        create_bp_collection_from_base_pairs_file(input_filename);

    // internal sequence dependence model
    if (m_cmd_line.DNA_model_type() == ForceFieldType::Internal)
        m_bp_collection.set_sequence_dependence_model(m_cmd_line.DNA_model());

    // external sequence dependence model
    else {

        // check for externnal force field
        DS_ASSERT(m_cmd_line.DNA_model_type() == ForceFieldType::External,
                  "fatal error in force field type");

        // filename
        const std::string ffield_filename = m_cmd_line.DNA_model();

        // binary archive
        InputArchive<CerealBinaryInput> bin_archive(ffield_filename);

        // loading
        StepParametersDB intrinsic_steps;
        ForceConstantsDB force_constants;
        bin_archive.load(intrinsic_steps);
        bin_archive.load(force_constants);

        // setup
        m_bp_collection.set_sequence_dependence_model(intrinsic_steps,
                                                      force_constants);

    };

    // frozen steps
    m_bp_collection.set_frozen_steps_domains(m_cmd_line.frozen_steps_list());

    // bp collection interface - free collection
    if (m_cmd_line.collection_type() == CollectionType::FreeCollection) {
        m_bp_intf_ptr = std::make_shared<FreeBpCollection>();
        m_bp_intf_ptr->set_bp_collection(m_bp_collection);
    }

    // bp collection interface - anchored collection
    else if (m_cmd_line.collection_type() ==
             CollectionType::AnchoredCollection) {
        m_bp_intf_ptr = std::make_shared<AnchoredBpCollection>();
        m_bp_intf_ptr->set_bp_collection(m_bp_collection);
    }

    // bp collection interface - clamped collection
    else if (m_cmd_line.collection_type() ==
             CollectionType::ClampedCollection) {
        m_bp_intf_ptr = std::make_shared<ClampedBpCollection>();
        m_bp_intf_ptr->set_bp_collection(m_bp_collection);
    }

    // bp collection interface - pulling collection
    else if (m_cmd_line.collection_type() ==
             CollectionType::PullCollection) {

        // pulling force
        Vector3 f;
        m_cmd_line.pulling_force(f);

        // temporary collection interface
        std::shared_ptr<PullingBpCollection> ptr =
        std::make_shared<PullingBpCollection>();
        ptr->set_bp_collection(m_bp_collection);
        ptr->set_pulling_force(f);

        // final interface
        m_bp_intf_ptr = std::move(ptr);

    }

    // checkpoint
    else {
        DS_ASSERT(false, "fatal error in collection type");
    };

    // electrostatics
    if (m_cmd_line.dh_electrostatics())
        m_bp_intf_ptr->toggle_on_electrostatics();
    else
        m_bp_intf_ptr->toggle_off_electrostatics();

};


// log file setup method
void emDNA_Application::setup_log_file(const std::string& output_name) {
    m_log_file.create_log_file(output_name + ".log");
};


// application setup message method
void emDNA_Application::application_setup_message() const {

    std::cout << "--- bp collection input\n";

    // input file info
    std::cout << "  bp collection created from input file: ";
    std::cout << m_cmd_line.input_file() << "\n";

    // size info
    std::cout << "  bp collection size: ";
    std::cout << m_bp_collection.n_of_base_pairs() << "-bp ";
    std::cout << "(" << m_bp_collection.n_of_bp_steps() << " steps)\n";

    // sequence dependence
    std::cout << "  DNA sequence-dependence model: ";
    if (m_cmd_line.DNA_model_type() == ForceFieldType::Internal)
        std::cout << m_cmd_line.DNA_model() << "\n";
    else {
        std::cout << "  to be loaded from file " + m_cmd_line.DNA_model();
        std::cout << "\n";
    };

    // frozen steps
    const std::vector<SizePair>& frozen_domains =
    m_bp_collection.frozen_steps_domains();
    std::cout << "  bp collection frozen steps: ";
    print_vector(frozen_domains);

    // electrostatics
    std::cout << "  DNA electrostatics included: ";
    if (m_cmd_line.dh_electrostatics())
        std::cout << "yes\n";
    else
        std::cout << "no\n";

    std::cout << "\n";

};


// pre-minimization message
void emDNA_Application::
minimization_pre_flight_message(const AlglibMinSettings_Ptr& mset) const {

    std::cout << "--- minimization setup\n";

    // minimization type
    std::cout << "  minimization type: ";
    std::cout << minimization_type_description(m_cmd_line.
                                               collection_type()) << "\n";

    // pulling force if needed
    if (m_cmd_line.collection_type() == CollectionType::PullCollection) {
        std::cout << "  pulling force: ";
        Vector3 f;
        DS_ASSERT(m_cmd_line.pulling_force(f),
                  "a pulling force is expected but not defined");
        std::cout << f << "\n";
    };

    // minim settings
    if (mset != nullptr) {
        std::cout << "  minimization settings overriden by command line:\n";
        std::cout << std::fixed << std::setprecision(REAL_WIDTH);
        std::cout << "    max-iterations = " << mset->_max_iterations << "\n";
        std::cout << "    threshold_dx = " << mset->_threshold_dx << "\n";
        std::cout << "    threshold_f = " << mset->_threshold_f << "\n";
        std::cout << "    threshold_g = " << mset->_threshold_g << "\n";
        std::cout << "    max-step-size = " << mset->_max_step_size << "\n";
    };

    std::cout << "\n";

};


// post-minimization message
void emDNA_Application::
minimization_post_flight_message(const MinimizationResults& res) const {

    std::cout << "--- minimization results\n";

    // initial and final value
    std::cout << "  energy initial value: ";
    std::cout << std::fixed << std::setprecision(REAL_WIDTH);
    std::cout << res._initial_E << "\n";
    std::cout << "  energy final value: ";
    std::cout << std::fixed << std::setprecision(REAL_WIDTH);
    std::cout << res._final_E << "\n";

    // gradient norm
    std::cout << "  gradient norm: ";
    std::cout << std::fixed << std::setprecision(REAL_WIDTH);
    std::cout << res._free_dofs_gradient.norm() << "\n";

    // number of iterations
    std::cout << "  # iterations: " << res._n_iterations << "\n";

    // return code
    std::cout << "  return code: " << res._return_code << "\n";

    std::cout << "\n";

};


// application tear down message
void emDNA_Application::apllication_tear_down_message() const {

};


// static method for minimization type description
std::string
emDNA_Application::minimization_type_description(const CollectionType&
                                                 coll_type) {

    std::string desc;

    // alglib free collection
    if (coll_type == CollectionType::FreeCollection)
        desc = "Alglib gradient-based minimization - free collection (DEBUG)";

    // alglib fixed origin
    if (coll_type == CollectionType::AnchoredCollection)
        desc = "Alglib gradient-based minimization - last bp origin fixed";

    // alglib fixed last bp
    if (coll_type == CollectionType::ClampedCollection)
        desc = "Alglib gradient-based minimization - last bp fixed";

    // alglib pulling
    if (coll_type == CollectionType::PullCollection)
        desc = "Alglib gradient-based minimization - pulling";

    return desc;

};


// ----


// entry point
int main(int argc, char* argv[]) {

    try {

        // new application
        emDNA_Application emDNA_app;

        // command line parsing
        emDNA_app.parse_command_line(argc, argv);

        // application setup
        emDNA_app.setup_application();

        // minimization
        emDNA_app.perform_minimization();

        // output
        emDNA_app.output_optimized_bp_collection();

    }
    
    catch(alglib::ap_error& e) {
        printf("error msg: %s\n", e.msg.c_str());
    }
    catch(DNASim_ExitException& e) {
        exit(e._exit_code);
    };

    return 0;
    
};

