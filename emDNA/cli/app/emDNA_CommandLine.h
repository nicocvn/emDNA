// emDNA_CommandLine function
// Nicolas Clauvelin


// command line parser based on TCLAP


#ifndef emDNA_emDNA_CommandLine_h
#define emDNA_emDNA_CommandLine_h


#include <TCLAP_Includes.h>
#include <emDNA_Includes.h>


// enum classes
enum class InputMode : Integer {
    BpList = 0,
    x3DNA_Bp = 1,
    x3DNA_Params = 2
};
enum class CollectionType : Integer {
    FreeCollection = 0,
    AnchoredCollection = 1,
    ClampedCollection = 2,
    PullCollection = 3
};
enum class ForceFieldType : Integer {
    Internal = 0,
    External = 1
};


class emDNA_CommandLine {


public:

    // constructors
    emDNA_CommandLine();
    ~emDNA_CommandLine();

    // command line parsing method
    void parse_command_line_input(int argc, char* argv[]);

    // input data report method
    void report_command_line_input_data(std::stringstream& output) const;

    // command line data accessors
    InputMode               input_mode() const;
    std::string             input_file() const;
    CollectionType          collection_type() const;
    bool                    pulling_force(Vector3& f) const;
    ForceFieldType          DNA_model_type() const;
    std::string             DNA_model() const;
    bool                    dh_electrostatics() const;
    std::string             output_name() const;
    std::vector<SizePair>   frozen_steps_list() const;
    bool                    minim_settings(AlglibMinSettings&
                                           minim_settings) const;
    bool                    quiet_flag() const;
    bool                    energy_progress_flag() const;
    Size                    output_progress_freq() const;


private:

    // command line arguments and options
    const std::vector<std::string> CommandLineParameters =
    std::vector<std::string>({
        "input_mode",
        "input_file",
        "collection_type",
        "pulling_force",
        "output_name",
        "DNA_model",
        "dh_electrostatics",
        "frozen_steps_list",
        "minim_settings",
        "quiet_flag",
        "energy_progress_flag",
        "output_progress_freq"
    });

    // private option setting method with assert checking
    template <class KeyType>
    void set_option_with_assert_check(const std::string& opt,
                                      KeyType opt_value) {
        DS_ASSERT(m_cl_options.set_option(opt, opt_value),
                  "setting non-existing command line option " + opt);
    };

    // private option finding method with assert checking
    template <class KeyType>
    void find_option_with_assert_check(const std::string& opt,
                                       KeyType& opt_value) const {
        DS_ASSERT(m_cl_options.get_option(opt, opt_value),
                  "requesting non-existing command line option " + opt);
    };

    // options setup private methods
    void setup_input_mode_and_file();
    void setup_collection_type();
    void setup_pulling_force();
    void setup_DNA_model();
    void setup_dh_electrostatics_flag();
    void setup_output_name();
    void setup_frozen_steps_list();
    void setup_minim_settings();
    void setup_quiet_flag();
    void setup_energy_progress_flag();
    void setup_output_progress_freq();

    // TCLAP command line object instantiation method
    // this is where the command line credits are setup
    static std::unique_ptr<TCLAP::CmdLine> emDNA_TCLAP_cmd_line_ptr();

    // command line arguments declaration method
    // this is where the arguments are defined and setup
    void declare_command_line_arguments(std::unique_ptr<TCLAP::CmdLine>&
                                        cmd_line_ptr);

    // options manager
    OptionsManager m_cl_options;

    // command line pointer object
    std::unique_ptr<TCLAP::CmdLine> m_cmd_ptr;


};


#endif  // emDNA_emDNA_CommandLine_h
