// CommandLineElements header
// Nicolas Clauvelin


#ifndef emDNA_CommandLineElements_h
#define emDNA_CommandLineElements_h


#include <TCLAP_Includes.h>


namespace {


    // input mode arguments declaration
    // one of those arguments is required
    TCLAP::ValueArg<std::string>
    x3DNAbp_input_arg("", "x3DNA-bp-input",
                      "x3DNA base pairs input file.",
                      true,
                      "",
                      "filename");
    TCLAP::ValueArg<std::string>
    x3DNAparams_input_arg("",
                          "x3DNA-bp-step-params-input",
                          "x3DNA bp step parameters input file.",
                          true,
                          "",
                          "filename");
    TCLAP::ValueArg<std::string>
    bplist_input_arg("", "bp-list-input",
                     "Base pairs input file (list format).",
                     true,
                     "",
                     "filename");

    // collection type switches declaration
    // one of those arguments is required
    TCLAP::SwitchArg
    free_collection_sw("", "free-collection",
                       "Minimization with no end conditions (for debug only).",
                       false);
    TCLAP::SwitchArg
    hold_last_origin_sw("", "hold-last-origin",
                        "Minimization with last origin imposed.",
                        false);
    TCLAP::SwitchArg
    hold_last_bp_sw("", "hold-last-bp",
                    "Minimization with last bp imposed.",
                    false);
    TCLAP::SwitchArg
    pull_last_bp_sw("", "pull-last-bp",
                    "Minimization with pulling force on last bp.",
                    false);

    // option pulling force
    TCLAP::ValueArg<std::string>
    pulling_force_arg("", "pulling-force",
                      "Pulling force.",
                      false,
                      "",// no default value
                      "{vector}");

    // option - sequence dependence
    TCLAP::ValueArg<std::string>
    DNA_model_arg("", "DNA-seqdep-model",
                  "DNA sequence-dependence model.",
                  false,
                  "IdealDNA",// default value
                  "string");
    TCLAP::ValueArg<std::string>
    external_DNA_model_arg("", "DNA-external-model",
                           "File containing a binary force field.",
                           false,
                           "",
                           "filename (binary)");

    // option - electrostatics
    TCLAP::SwitchArg
    add_dh_electrostatics_sw("", "dh-electrostatics",
                             "Add electrostatics interactions to the objective "
                             "function.",
                             false);

    // option - output name
    TCLAP::ValueArg<std::string>
    output_name_arg("", "output-name",
                    "Output prefix name. Default is 'emDNA_minim'.",
                    false,
                    "emDNA_minim",// default value
                    "string");

    // option - frozen steps
    TCLAP::ValueArg<std::string>
    frozen_steps_arg("", "frozen-steps",
                     "List of ranges of frozen steps n1:n2 (numbering starts "
                     "at 1 and the range is inclusive).",
                     false,
                     "",
                     "{list}");

    // option - minimization settings
    TCLAP::ValueArg<std::string>
    minim_settings_arg("", "minim-settings",
                       "{max-iterations, dx, f, g, max-step-size}.",
                       false,
                       "",
                       "{vector}");

    // option - quiet output
    TCLAP::SwitchArg
    quiet_sw("", "quiet",
             "Suppress console output.",
             false);

    // option - energy progress
    TCLAP::SwitchArg
    energy_progress_sw("", "energy-progress",
                       "Display the current energy during the minimization.",
                       false);

    // option - output progress
    TCLAP::ValueArg<std::string>
    output_progress_arg("",
                        "output-progress",
                        "Frequency for configurations output "
                        "(output file is tmp_confs.txt). "
                        "Zero frequencydisables output.",
                        false,
                        "0",
                        "integer");


    // list of all ValueArg and Switch Arg
    std::vector<TCLAP::ValueArg<std::string>*> ValueArg_list({
        &x3DNAbp_input_arg, &x3DNAparams_input_arg, &bplist_input_arg,
        &pulling_force_arg,
        &DNA_model_arg,
        &external_DNA_model_arg,
        &output_name_arg,
        &frozen_steps_arg,
        &minim_settings_arg,
        &output_progress_arg
    });
    std::vector<TCLAP::SwitchArg*> SwitchArg_list({
        &free_collection_sw, &hold_last_origin_sw, &hold_last_bp_sw,
            &pull_last_bp_sw,
        &add_dh_electrostatics_sw,
        &quiet_sw,
        &energy_progress_sw
    });


}


#endif  //emDNA_CommandLineElements_h
