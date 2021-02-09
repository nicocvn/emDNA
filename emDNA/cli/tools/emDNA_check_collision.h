// emDNA_check_collision tool
// Nicolas Clauvelin


#ifndef emDNA_emDNA_check_collision_h
#define emDNA_emDNA_check_collision_h


#include <emDNA_Includes.h>


// data structure
struct ParserData {
    std::string _filename;
    bool _input_x3DNAbp = false;
    bool _input_x3DNAparams = false;
    bool _input_bplist = false;
    bool _closed = false;
    Real _DNA_radius;
};


// entry point
Integer main(Integer argc, char* argv[]);


// command line parsing option
ParserData parse_command_line(Integer argc, char* argv[]);


// collision checking methods
bool linear_fragment_collision_check(const BpCollection& bp_coll,
                                     const Real& DNA_radius);
bool closed_fragment_collision_check(const BpCollection& bp_coll,
                                     const Real& DNA_radius);


#endif  // emDNA_emDNA_check_collision_h
