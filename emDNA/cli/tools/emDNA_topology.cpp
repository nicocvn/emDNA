// emDNA_topology tool
// Nicolas Clauvelin


#include <TCLAP_Includes.h>
#include <emDNA.h>
#include <emDNA_topology.h>


// entry point
Integer main(Integer argc, char* argv[]) {

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
        }
        else if (parser_data._input_x3DNAparams) {
            bp_collection =
            BpCollectionFactory::
            create_bp_collection_from_x3DNA_bp_step_params_file(parser_data.
                                                                _filename);
        }
        else {
            bp_collection =
            BpCollectionFactory::
            create_bp_collection_from_base_pairs_file(parser_data._filename);
        };

        // virtual last bp
        // if set we remove the last base pair
        if (parser_data._virtual_last_bp) {
            std::vector<BasePair> all_base_pairs = bp_collection.base_pairs();
            all_base_pairs.pop_back();
            bp_collection =
            BpCollection::collection_from_base_pairs(all_base_pairs);
        };

        // bp collection statistics
        std::cout << "--- bp collection input ---\n";
        std::cout << "  bp collection created from input file: ";
        std::cout << parser_data._filename << "\n";
        std::cout << "  bp collection size: ";
        std::cout << bp_collection.n_of_base_pairs() << "-bp ";
        std::cout << "(" << bp_collection.n_of_bp_steps() << " steps)\n";
        if (parser_data._virtual_last_bp)
            std::cout << "  last bp dropped because it is virtual\n";
        std::cout << "\n";

        // writhe
        const Real bp_coll_Wr = bp_collection_writhe(bp_collection);

        // total twist
        const std::vector<Real> bp_coll_Tw =
        bp_collection_total_twist_density(bp_collection);

        // direct Lk
        const Real bp_coll_Lk = bp_collection_direct_lk(bp_collection);

        // twist density
        if (parser_data._twist_density) {
            std::cout << "--- twist density per step ---\n";
            for (Size i=0; i<bp_coll_Tw.size(); ++i) {
                std::cout << i << "\t";
                std::cout << std::setprecision(REAL_WIDTH) << bp_coll_Tw[i]
                    << "\n";
            };
            std::cout << "\n";
        };

        // total twist
        Real Tw = FLOAT_INIT;
        std::for_each(bp_coll_Tw.begin(), bp_coll_Tw.end(), [&Tw](const Real&
                                                                  dtw){
            Tw += dtw;
        });
        Tw /= (2*F_PI);

        // topology results
        std::cout << "--- toplogy results ---\n";
        std::cout << "Wr = " << std::setprecision(REAL_WIDTH) << bp_coll_Wr
            << "\n";
        std::cout << "Tw = " << std::setprecision(REAL_WIDTH) << Tw << "\n";
        std::cout << "Lk = " << std::setprecision(REAL_WIDTH) << bp_coll_Lk
            << "\n";
        std::cout << "Lk_eps = " << std::setprecision(REAL_WIDTH) <<
            bp_coll_Lk-(bp_coll_Wr+Tw);
        std::cout << "\n";

        return 0;

    }

    catch(DNASim_ExitException& e) {
        exit(e._exit_code);
    };

};


// command line parsing option
ParserData parse_command_line(Integer argc, char* argv[]) {

    // command line parser
    TCLAP::CmdLine cmd_line("emDNA_check_collision - Nicolas Clauvelin, "
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
                                              "bp list input file.",
                                              true,
                                              "",
                                              "string");
    std::vector<TCLAP::Arg*> xor_input;
    xor_input.push_back(&input_x3DNAbp);
    xor_input.push_back(&input_x3DNAparams);
    xor_input.push_back(&input_bplist);
    cmd_line.xorAdd(xor_input);

    // flags
    TCLAP::SwitchArg last_is_virtual("", "virtual-last-bp",
                                     "Assume the last bp is virtual "
                                     "(useful for minicircles).",
                                     false);
    cmd_line.add(last_is_virtual);
    TCLAP::SwitchArg twist_density("", "twist-density",
                                   "Output twist density results.",
                                   false);
    cmd_line.add(twist_density);

    // parsing
    cmd_line.parse(argc, argv);

    // input settings
    ParserData parser_data;
    if (input_x3DNAbp.isSet()) {
        parser_data._input_x3DNAbp = true;
        parser_data._filename = input_x3DNAbp.getValue();
    }
    else if (input_x3DNAparams.isSet()){
        parser_data._input_x3DNAparams = true;
        parser_data._filename = input_x3DNAparams.getValue();
    }
    else {
        parser_data._input_bplist = true;
        parser_data._filename = input_bplist.getValue();
    };

    // flags values
    if (last_is_virtual.isSet()) {
        parser_data._virtual_last_bp = true;
    };
    if (twist_density.isSet()) {
        parser_data._twist_density = true;
    };

    return parser_data;
    
};


// bp collection writhing number
// always assumes a closed collection
// remove any virtual point before hand
Real bp_collection_writhe(const BpCollection& bp_coll) {

    // collection size
    const Size n_bp = bp_coll.n_of_base_pairs();

    // collection centerline
    CurvePoints centerline;
    centerline.reserve(n_bp);
    for (Size i=0; i<n_bp; ++i)
        centerline.push_back(bp_coll.base_pair(i).origin());

    return CurveTopology::curve_writhing_number(centerline);

};


// bp collection direct linking number
// always assumes a closed collection
// remove any virtual point before hand
Real bp_collection_direct_lk(const BpCollection& bp_coll) {

    // collection size
    const Size n_bp = bp_coll.n_of_base_pairs();

    // collection centerline
    CurvePoints centerline;
    centerline.reserve(n_bp);
    for (Size i=0; i<n_bp; ++i)
        centerline.push_back(bp_coll.base_pair(i).origin());

    // collection edge curve
    // we pick an arbitrary small value for the linking number to avoid
    // collision issues
    CurvePoints edge;
    edge.reserve(n_bp);
    for (Size i=0; i<n_bp; ++i)
#define EPS_VEC Real(0.01)
        edge.push_back(bp_coll.base_pair(i).origin()+
                       EPS_VEC*bp_coll.base_pair(i).axis(I));
#undef EPS_VEC

    // discrete ribbon
    DiscreteRibbon ribb(centerline, edge);

    // Lk
    return ribb.ribbon_linking_number();
    
};


// bp collection total twist
// the modification for the first and last steps are performed in BpTwistDensity
std::vector<Real> bp_collection_total_twist_density(const BpCollection&
                                                    bp_coll) {

    // cylization
    // we add the first two base pairs at the end
    // we add the last base pair at the beginning
    std::vector<BasePair> base_pairs = bp_coll.base_pairs();

    // twist results
    std::vector<Vector3> twist_res =
    BpTwistDensity::compute_twist_density(base_pairs);

    // size checking
    DS_ASSERT(twist_res.size() == bp_coll.n_of_base_pairs(),
              "sizing problem in twist density computation");

    // total twist
    // we discard the first and last values
    std::vector<Real> Tw;
    for (Size i=0; i<twist_res.size(); ++i)
            Tw.push_back(twist_res[i][Z]);

    return Tw;

};
