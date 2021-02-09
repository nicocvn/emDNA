// emDNA_CommandLine function
// Nicolas Clauvelin


#include <emDNA.h>
#include <CommandLineElements.h>
#include <emDNA_CommandLine.h>


// boolean strings
#define TRUE_STR "1"
#define FALSE_STR "0"


// override usage output method to display list of sequence dependence models
class emDNAOutput : public TCLAP::StdOutput {
public:
    virtual void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e) {
        TCLAP::StdOutput::failure(c, e);
    };

    virtual void usage(TCLAP::CmdLineInterface& c) {

        TCLAP::StdOutput::usage(c);

        // force field inventory
        std::vector<std::string> inventory;
        for (auto it = SequenceDependenceModelList.begin();
             it !=SequenceDependenceModelList.end(); ++it)
            inventory.push_back(it->first
                                + "\n\t\t"
                                + it->second._description);

        // output
        std::cout << "Implemented force fields:\n";
        for (auto s : inventory)
            std::cout << "\t" << s << "\n";
        std::cout << "\n";

    };

    virtual void version(TCLAP::CmdLineInterface& c) {
        TCLAP::StdOutput::version(c);
    };

};


// class default constructor
emDNA_CommandLine::emDNA_CommandLine() :
m_cl_options(),
m_cmd_ptr(nullptr) {

    // set options list
    m_cl_options.set_options_list(CommandLineParameters);

    // command line object
    m_cmd_ptr = emDNA_TCLAP_cmd_line_ptr();


    // custom output
    m_cmd_ptr->setOutput(new emDNAOutput());

};


// class destructor
emDNA_CommandLine::~emDNA_CommandLine() {
    pointer_delete(m_cmd_ptr->getOutput());
};


// command line parsing method
void emDNA_CommandLine::parse_command_line_input(int argc, char* argv[]) {

    // command line arguments declaration
    declare_command_line_arguments(m_cmd_ptr);

    // command line parsing
    m_cmd_ptr->parse(argc, argv);

    // options setup
    setup_input_mode_and_file();
    setup_collection_type();
    setup_pulling_force();
    setup_DNA_model();
    setup_output_name();
    setup_frozen_steps_list();
    setup_minim_settings();
    setup_quiet_flag();
    setup_energy_progress_flag();
    setup_output_progress_freq();
    setup_dh_electrostatics_flag();

    // consistentcy checks
    Vector3 f;
    if (!pull_last_bp_sw.isSet() && pulling_force(f)) {
        std::cout << "a pulling force can only be set in pull-last-bp mode\n";
        std::cout << "./emDNA -h for help\n";
        std::cout << "\n\n";
        throw DNASim_ExitException(-1);
    };
    if (pull_last_bp_sw.isSet() && !pulling_force(f)) {
        std::cout << "pull-last-bp mode without pulling force\n";
        std::cout << "./emDNA -h for help\n";
        std::cout << "\n\n";
        throw DNASim_ExitException(-1);
    };

};


// input data report method
void emDNA_CommandLine::
report_command_line_input_data(std::stringstream& output) const {

    output << std::boolalpha;

    // loop over all arguments
    for (auto it = ValueArg_list.begin(); it != ValueArg_list.end(); ++it) {
        output << (*it)->longID() << ": ";
        output << (*it)->getValue() << "\n";
    };
    for (auto it = SwitchArg_list.begin(); it != SwitchArg_list.end(); ++it) {
        output << (*it)->longID() << ": ";
        output << (*it)->getValue() << "\n";
    };

};


// command line data accessors
InputMode emDNA_CommandLine::input_mode() const {

    // get enum integer value
    Integer input_mode_val;
    find_option_with_assert_check("input_mode", input_mode_val);

    // cast to enum type
    return static_cast<InputMode>(input_mode_val);

};
std::string emDNA_CommandLine::input_file() const {

    // get input file string
    std::string input_file;
    find_option_with_assert_check("input_file", input_file);

    return input_file;

};
CollectionType emDNA_CommandLine::collection_type() const {

    // get enum integer value
    Integer collection_type_val;
    find_option_with_assert_check("collection_type", collection_type_val);

    // cast to enum type
    return static_cast<CollectionType>(collection_type_val);

};
bool emDNA_CommandLine::pulling_force(Vector3& f) const {

    // vector string
    std::string f_str;
    find_option_with_assert_check("pulling_force", f_str);

    // if string is zero length there is no settings
    if (f_str.length() == 0)
        return false;

    // otherwise parsing
    else {
        f = Vector3(f_str);
        return true;
    };

};
ForceFieldType emDNA_CommandLine::DNA_model_type() const {

    // get model string
    std::string model_str;
    find_option_with_assert_check("DNA_model", model_str);

    // check if it contains 'external' keyword
    Size found = model_str.find("external");

    if (found != std::string::npos)
        return ForceFieldType::External;
    else
        return ForceFieldType::Internal;

};
std::string emDNA_CommandLine::DNA_model() const {

    // get model string
    std::string model_str;
    find_option_with_assert_check("DNA_model", model_str);

    // check for external keyword
    Size found = model_str.find("external");
    if (found != std::string::npos) {
        std::vector<std::string> tokens =
        EnhancedString::tokenize_string(model_str, ':');
        return tokens.back();
    };

    return model_str;

};
bool emDNA_CommandLine::dh_electrostatics() const {

    // get boolean value
    bool dh_electrostatics_flag;
    find_option_with_assert_check("dh_electrostatics", dh_electrostatics_flag);

    return dh_electrostatics_flag;

};
std::string emDNA_CommandLine::output_name() const {

    // get model string
    std::string output_name_str;
    find_option_with_assert_check("output_name", output_name_str);

    return output_name_str;
    
};
std::vector<SizePair> emDNA_CommandLine::frozen_steps_list() const {

    // get frozen steps list string
    std::string frozen_steps_str;
    find_option_with_assert_check("frozen_steps_list", frozen_steps_str);

    // parsing
    return emDNA_Utils::parse_frozen_steps_list(frozen_steps_str);

};
bool emDNA_CommandLine::minim_settings(AlglibMinSettings& minim_settings) const
{

    // minim settings string
    std::string settings_str;
    find_option_with_assert_check("minim_settings", settings_str);

    // if string is zero length there is no settings
    if (settings_str.length() == 0)
        return false;

    // otherwise parsing
    else {
        minim_settings = emDNA_Utils::parse_minim_settings(settings_str);
        return true;
    };

};
bool emDNA_CommandLine::quiet_flag() const {

    // get boolean value
    bool quiet_flag_bool;
    find_option_with_assert_check("quiet_flag", quiet_flag_bool);

    return quiet_flag_bool;

};
bool emDNA_CommandLine::energy_progress_flag() const {

    // get boolean value
    bool energy_progress_flag_bool;
    find_option_with_assert_check("energy_progress_flag",
                                  energy_progress_flag_bool);

    return energy_progress_flag_bool;
    
};
Size emDNA_CommandLine::output_progress_freq() const {

    // get frequency value
    Size freq;
    find_option_with_assert_check("output_progress_freq", freq);

    return freq;

};


// options setup private method - input mode
void emDNA_CommandLine::setup_input_mode_and_file() {

    // check each argument and adjust options accordingly
    if (x3DNAbp_input_arg.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)InputMode::x3DNA_Bp);
        set_option_with_assert_check("input_mode", s);
        set_option_with_assert_check("input_file",
                                     x3DNAbp_input_arg.getValue());
    }
    else if (x3DNAparams_input_arg.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)InputMode::x3DNA_Params);
        set_option_with_assert_check("input_mode", s);
        set_option_with_assert_check("input_file",
                                     x3DNAparams_input_arg.getValue());
    }
    else if (bplist_input_arg.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)InputMode::BpList);
        set_option_with_assert_check("input_mode", s);
        set_option_with_assert_check("input_file",
                                     bplist_input_arg.getValue());
    }
    else {
        DS_ASSERT(false,
                  "fatal error in input mode parsing");
    };

};


// options setup private method - collection type
void emDNA_CommandLine::setup_collection_type() {

    // check each argument and adjust option accordingly
    if (free_collection_sw.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)
                                          CollectionType::FreeCollection);
        set_option_with_assert_check("collection_type", s);
    }
    else if (hold_last_origin_sw.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)
                                          CollectionType::AnchoredCollection);
        set_option_with_assert_check("collection_type", s);
    }
    else if (hold_last_bp_sw.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)
                                          CollectionType::ClampedCollection);
        set_option_with_assert_check("collection_type", s);
    }
    else if (pull_last_bp_sw.isSet()) {
        const std::string s =
        EnhancedString::convert_to_string((Integer)
                                          CollectionType::PullCollection);
        set_option_with_assert_check("collection_type", s);
    }
    else {
        DS_ASSERT(false,
                  "fatal error in collection type parsing");
    };
    
};


// options setup private method - pulling force
void emDNA_CommandLine::setup_pulling_force() {
    set_option_with_assert_check("pulling_force",
                                 pulling_force_arg.getValue());
};


// options setup private method - DNA model
void emDNA_CommandLine::setup_DNA_model() {
    if (DNA_model_arg.isSet())
        set_option_with_assert_check("DNA_model", DNA_model_arg.getValue());
    else if (external_DNA_model_arg.isSet())
        set_option_with_assert_check("DNA_model",
                                     "external:" +
                                     external_DNA_model_arg.getValue());
    else {
        DS_ASSERT(false,
                  "fatal error in DNA model parsing");
    }
};


// options setup private method - dh electrostatics flag
void emDNA_CommandLine::setup_dh_electrostatics_flag() {
    if (add_dh_electrostatics_sw.isSet())
        set_option_with_assert_check("dh_electrostatics", TRUE_STR);
    else
        set_option_with_assert_check("dh_electrostatics", FALSE_STR);
};


// options setup private method - output model
void emDNA_CommandLine::setup_output_name() {
    set_option_with_assert_check("output_name", output_name_arg.getValue());
};


// options setup private method - frozen steps list
void emDNA_CommandLine::setup_frozen_steps_list() {
    set_option_with_assert_check("frozen_steps_list",
                                 frozen_steps_arg.getValue());
};


// options setup private method - minim settings
void emDNA_CommandLine::setup_minim_settings() {
    set_option_with_assert_check("minim_settings",
                                 minim_settings_arg.getValue());
};


// options setup private method - quiet flag
void emDNA_CommandLine::setup_quiet_flag() {
    if (quiet_sw.isSet())
        set_option_with_assert_check("quiet_flag", TRUE_STR);
    else
        set_option_with_assert_check("quiet_flag", FALSE_STR);
};


// options setup private method - energy progress flag
void emDNA_CommandLine::setup_energy_progress_flag() {
    if (energy_progress_sw.isSet())
        set_option_with_assert_check("energy_progress_flag", TRUE_STR);
    else
        set_option_with_assert_check("energy_progress_flag", FALSE_STR);
};


// options setup private method - output progress freq
void emDNA_CommandLine::setup_output_progress_freq() {
    set_option_with_assert_check("output_progress_freq",
                                 EnhancedString::
                                 convert_to_string(output_progress_arg.
                                                   getValue()));
};


// TCLAP command line object instantiation method
// this is where the command line are credits are setup
std::unique_ptr<TCLAP::CmdLine> emDNA_CommandLine::emDNA_TCLAP_cmd_line_ptr() {

    // new pointer
#define emDNA_VERSION "0.x"
    std::unique_ptr<TCLAP::CmdLine>
    cmd_line_ptr(new TCLAP::CmdLine("emDNA - Nicolas Clauvelin, "
                                    "Rutgers University",
                                    '=',
                                    emDNA_VERSION));
#undef emDNA_VERSION

    return cmd_line_ptr;

};


// command line arguments declaration method
// this is where the arguments are defined and setup
void emDNA_CommandLine::
declare_command_line_arguments(std::unique_ptr<TCLAP::CmdLine>& cmd_line_ptr) {

    // input arguments declaration
    // one of those arguments is required
    std::vector<TCLAP::Arg*> xor_input;
    xor_input.push_back(&x3DNAbp_input_arg);
    xor_input.push_back(&x3DNAparams_input_arg);
    xor_input.push_back(&bplist_input_arg);
    cmd_line_ptr->xorAdd(xor_input);

    // collection type switches declaration
    // one of those arguments is required
    std::vector<TCLAP::Arg*> xor_switches;
    xor_switches.push_back(&free_collection_sw);
    xor_switches.push_back(&hold_last_origin_sw);
    xor_switches.push_back(&hold_last_bp_sw);
    xor_switches.push_back(&pull_last_bp_sw);
    cmd_line_ptr->xorAdd(xor_switches);

    // pulling force
    cmd_line_ptr->add(&pulling_force_arg);

    // option DNA model and external DNA models
    std::vector<TCLAP::Arg*> xor_DNA_model;
    xor_DNA_model.push_back(&DNA_model_arg);
    xor_DNA_model.push_back(&external_DNA_model_arg);
    cmd_line_ptr->xorAdd(xor_DNA_model);

    // option output name
    cmd_line_ptr->add(&output_name_arg);

    // option frozen steps
    cmd_line_ptr->add(&frozen_steps_arg);

    // option minim settings
    cmd_line_ptr->add(&minim_settings_arg);

    // option quiet output flag
    cmd_line_ptr->add(&quiet_sw);

    // option energy progress
    cmd_line_ptr->add(&energy_progress_sw);

    // option output progress
    cmd_line_ptr->add(&output_progress_arg);

    // option dh electrostatics flag
    cmd_line_ptr->add(&add_dh_electrostatics_sw);

};


#undef TRUE_STR
#undef FALSE_STR
