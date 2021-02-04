// emDNA_lego_gen tool
// Nicolas Clauvelin


#include <TCLAP_Includes.h>
#include <emDNA.h>
#include <emDNA_lego_gen.h>


// internal data for lego type translation
const std::map<std::string,LegoType> LegoTypeTranslator = {
    {"Bender", LegoType::Bender},
    {"Twister", LegoType::Twister},
    {"TwistedBender", LegoType::TwistedBender}
};


// entry point
int main(int argc, char* argv[]) {

    try {

        // command line parsing
        CLData cl_data = parse_command_line(argc, argv);

        // LEGO protein generation
        auto lego_ptr = build_lego_protein(cl_data);

        // output
        const std::vector<BasePair>& bps = lego_ptr->base_pairs();
        for (const BasePair& bp : bps)
            std::cout << bp << "\n";

    }

    catch (DNASim_ExitException& e) {
        exit(e._exit_code);
    };

};


// command line parsing option
CLData parse_command_line(int argc, char* argv[]) {

    // command line parser
    TCLAP::CmdLine cmd_line("emDNA_lego_gen - Nicolas Clauvelin, "
                            "Rutgers University",
                            '=',
                            "");

    // input arguments
    TCLAP::ValueArg<std::string> lego_type("", "lego-type",
                                           "LEGO protein type",
                                           true,
                                           "",
                                           "string");
    TCLAP::ValueArg<std::string> lego_size("", "lego-bp-size",
                                           "Size in bp of the LEGO protein.",
                                           true,
                                           "",
                                           "string");
    TCLAP::ValueArg<std::string> bending_angle("", "bending-angle",
                                               "Bending angle for LEGO "
                                               "protein.",
                                               false,
                                               "",
                                               "string");
    TCLAP::ValueArg<std::string> twisting_angle("", "twisting-angle",
                                                "Twisting angle for LEGO "
                                                "protein.",
                                                false,
                                                "",
                                                "string");

    // parsing
    cmd_line.add(lego_type);
    cmd_line.add(lego_size);
    cmd_line.add(bending_angle);
    cmd_line.add(twisting_angle);
    cmd_line.parse(argc, argv);

    // consistency check
    auto it = LegoTypeTranslator.find(lego_type.getValue());
    DS_ASSERT(it != LegoTypeTranslator.end(),
              "wrong LEGO protein type; "
              "types = {Bender, Twister, TwistedBender}");
    if (it->second == LegoType::Bender)
        DS_ASSERT(bending_angle.isSet(),
                  "bending angle required");
    if (it->second == LegoType::Twister)
        DS_ASSERT(twisting_angle.isSet(),
                  "twisting angle required");
    if (it->second == LegoType::TwistedBender)
        DS_ASSERT(twisting_angle.isSet() && bending_angle.isSet(),
                  "bending and twisting angles required");

    // input settings
    CLData cl_data;
    cl_data._lego_type = it->second;
    cl_data._n_bp = EnhancedString::convert_from_string<Size>(lego_size.
                                                              getValue());
    cl_data._bending_angle =
    EnhancedString::convert_from_string<Real>(bending_angle.getValue());
    cl_data._added_twist =
    EnhancedString::convert_from_string<Real>(twisting_angle.getValue());

    return cl_data;

};


// LEGO protein factory
std::unique_ptr<LegoProtein> build_lego_protein(const CLData& cl_data) {

    if (cl_data._lego_type == LegoType::Bender) {
        return
        std::unique_ptr<LegoBender>(new LegoBender(cl_data._n_bp,
                                                   cl_data._bending_angle));
    }
    else if (cl_data._lego_type == LegoType::Twister) {
        return
        std::unique_ptr<LegoTwister>(new LegoTwister(cl_data._n_bp,
                                                     cl_data._added_twist));
    }
    else if (cl_data._lego_type == LegoType::TwistedBender) {
        return
        std::unique_ptr<LegoTwistedBender>
        (new LegoTwistedBender(cl_data._n_bp,
                               cl_data._added_twist,
                               cl_data._bending_angle));
    }
    else {
        DS_ASSERT(false, "LEGO factory fatal error");
    };

};
