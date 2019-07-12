// emDNA_lego_gen tool
// Nicolas Clauvelin


#ifndef emDNA_emDNA_lego_gen_h
#define emDNA_emDNA_lego_gen_h


#include <emDNA_Includes.h>
class LegoProtein;


// lego type enum
enum class LegoType : Integer {
    Bender = 0,
    Twister = 1,
    TwistedBender = 2
};

// data structure
struct CLData {
    LegoType _lego_type;
    Size _n_bp;
    Real _bending_angle;
    Real _added_twist;
};

// entry point
int main(int argc, char* argv[]);

// command line parsing option
CLData parse_command_line(int argc, char* argv[]);

// LEGO protein factory
std::unique_ptr<LegoProtein> build_lego_protein(const CLData& cl_data);


#endif  // emDNA_emDNA_lego_gen_h
