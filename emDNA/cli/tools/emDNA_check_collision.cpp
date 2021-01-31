// emDNA_check_collision tool
// Nicolas Clauvelin


#include <TCLAP_Includes.h>
#include <BpCollection.h>
#include <BpCollectionFactory.h>
#include <emDNA_check_collision.h>


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

        // DNA radius
        const Real DNA_radius = parser_data._DNA_radius;

        // collision checking
        bool collision;
        if (parser_data._closed)
            collision = closed_fragment_collision_check(bp_collection,
                                                        DNA_radius);
        else
            collision = linear_fragment_collision_check(bp_collection,
                                                        DNA_radius);

        if (collision) {
            exit(1);
        }
        else {
            exit(0);
        };

    }

    catch (DNASim_ExitException& e) {
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

    // option - output name
    TCLAP::ValueArg<std::string> dna_radius("", "DNA-radius",
                                                "DNA radius for collision "
                                                "checking.",
                                                false,
                                                "10.0",    // default value
                                                "string");
    cmd_line.add(dna_radius);

    // option - closed fragment
    TCLAP::SwitchArg closed_frag("", "closed",
                                 "cyclic collection.");
    cmd_line.add(closed_frag);

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

    // DNA radius
    parser_data._DNA_radius =
    EnhancedString::convert_from_string<Real>(dna_radius.getValue());

    // cyclic collection
    if (closed_frag.isSet())
        parser_data._closed = true;

    return parser_data;

};


// collision checking methods
bool linear_fragment_collision_check(const BpCollection& bp_coll,
                                     const Real& DNA_radius) {

    // base pairs
    const std::vector<BasePair>& bps = bp_coll.base_pairs();

    // point cloud kd tree for all base pairs
    PointCloudKdTree<BasePair> bp_cloud_kdtree(bps);

    // query for each base pair
    const Size n_bps = bp_coll.n_of_base_pairs();
    for (Size bp_idx = 0; bp_idx<n_bps; ++bp_idx) {

        // radius search for points within twice the DNA radius from the bp
        // origin
        // the search is not sorted
        // we estimate at 10 the number of neighbors ... why not?
        RadiusResults radius_res =
        bp_cloud_kdtree.radius_search_from_point(bps[bp_idx].origin(),
                                                 2*DNA_radius,
                                                 10,
                                                 false);

        // results filter
        // for each found neighbor we check the sequential distance and reject
        // if the indices are less than some threshold apart
        const Size n_res = radius_res.size();
        for (Size i=0; i<n_res; ++i) {
            const Integer seq_dist =
            std::fabs(static_cast<Integer>(bp_idx)-
                      static_cast<Integer>(radius_res[i].first));
#define SEQ_THRESHOLD 10
            if (seq_dist < SEQ_THRESHOLD)
#undef SEQ_THRESHOLD
                continue;       // sequentially close neighbors
            else
                return true;    // real collision
        };

    };

    // no collision detected
    return false;

};


bool closed_fragment_collision_check(const BpCollection& bp_coll,
                                     const Real& DNA_radius) {

    // base pairs
    const std::vector<BasePair>& bps = bp_coll.base_pairs();

    // point cloud kd tree for all base pairs
    PointCloudKdTree<BasePair> bp_cloud_kdtree(bps);

    // query for each base pair
    const Size n_bps = bp_coll.n_of_base_pairs();
    for (Size bp_idx = 0; bp_idx<n_bps; ++bp_idx) {

        // radius search for points within twice the DNA radius from the bp
        // origin
        // the search is not sorted
        // we estimate at 10 the number of neighbors ... why not?
        RadiusResults radius_res =
        bp_cloud_kdtree.radius_search_from_point(bps[bp_idx].origin(),
                                                 2*DNA_radius,
                                                 10,
                                                 false);

        // results filter
        // for each found neighbor we check the sequential distance and reject
        // if the indices are less than some threshold apart
        const Size n_res = radius_res.size();
        for (Size i=0; i<n_res; ++i) {
            const Integer seq_dist =
            std::fabs(static_cast<Integer>(bp_idx)-
                      static_cast<Integer>(radius_res[i].first));
            const Integer comp_dist = n_bps-seq_dist;
            const Integer min_dist = std::min(seq_dist, comp_dist);
#define SEQ_THRESHOLD 10
            if (min_dist < SEQ_THRESHOLD)
#undef SEQ_THRESHOLD
                continue;       // sequentially close neighbors
            else {
                std::cout << bp_idx << "\t" << radius_res[i].first << "\n";
                return true;    // real collision
            }
        };

    };

    // no collision detected
    return false;

};
