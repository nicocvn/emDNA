// emDNA_topology tool
// Nicolas Clauvelin


// remarks on topology computation
// - the writing number of a bp collection is performed under the assumption
//  that the collection forms a closed centerline (be careful to remove any
//  virtual points as it will cause divergences),
// - same remark as for the writhing number holds for the direct linking number,
// - same remark as for the writhing number holds for the twist related results.

// about twist density:
// the collection is modified to make it possible to get twist density results
// for the first and last steps; if those values are not needed just discard
// them; keep in my mind that these values are used to compute the total twist
// ...
// so for an open collection one should output the twist density per step and
// make the summation in which the first and last step are thrown away
// ...
// for an open colleciton it is also possible to add additional bp to get values


#ifndef emDNA_emDNA_topology_h
#define emDNA_emDNA_topology_h


#include <emDNA_Includes.h>


// data structure
struct ParserData {
    std::string _filename;
    bool _input_x3DNAbp = false;
    bool _input_x3DNAparams = false;
    bool _input_bplist = false;
    bool _virtual_last_bp = false;
    bool _twist_density = false;
};


// entry point
Integer main(Integer argc, char* argv[]);


// command line parsing option
ParserData parse_command_line(Integer argc, char* argv[]);


// topology methods
Real bp_collection_writhe(const BpCollection& bp_coll);
Real bp_collection_direct_lk(const BpCollection& bp_coll);
std::vector<Real> bp_collection_total_twist_density(const BpCollection&
                                                    bp_coll);


#endif  // emDNA_emDNA_topology_h
