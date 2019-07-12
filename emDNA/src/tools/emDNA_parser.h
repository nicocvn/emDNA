// emDNA_parser tool
// Nicolas Clauvelin


#ifndef emDNA_emDNA_parser_h
#define emDNA_emDNA_parser_h


#include <emDNA_Includes.h>


// data structure
struct ParserData {
    std::string _filename;
    bool _input_x3DNAbp;
    bool _input_x3DNAparams;
    bool _input_bp_list;
    bool _output_bp_list;
    bool _output_x3DNAbp;
    bool _output_x3DNAparams;
    std::string _bp_list_sequence;
};


// entry point
int main(int argc, char* argv[]);


// command line parsing option
ParserData parse_command_line(int argc, char* argv[]);


#endif  // emDNA_emDNA_parser_h
