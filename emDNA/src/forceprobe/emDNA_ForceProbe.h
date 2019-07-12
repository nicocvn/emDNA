// emDNA_ForceProbe app
// Nicolas Clauvelin


// emDNA_ForceProbe input file example
//
// format is [ bplist | x3DNAbp | x3DNAparams ]
//
// base_collection = format:filename
// base_bound_domains = {n1_a:n1_b, n2_a:n2_b, ...} [{}]
// base_seqdep_model = model_name [IdealDNA]
// force_probe = {force_ini, force_fin, n_points}
// force_probe_direction = {X,Y,Z}


#ifndef emDNA_emDNA_ForceProbe_h
#define emDNA_emDNA_ForceProbe_h


#include <emDNA_Includes.h>
class BpCollection;


// entry point
int main(int argc, char* argv[]);

// input parsing function
struct AppData {
    std::string _base_input;
    std::string _base_bound_domains;
    std::string _base_seqdep_model;
    std::string _force_probe;
    std::string _force_probe_dir;
};
AppData parse_input(int argc, char* argv[]);

// bp collection creation function
BpCollection create_bp_collection(const std::string& input_string);

// inspector message function
void print_inspector_output(const BpCollection& base_collection,
                            const ForceRamp& force_ramp);


#endif  // emDNA_emDNA_ForceProbe_h

