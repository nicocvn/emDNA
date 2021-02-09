// emProBind app
// Nicolas Clauvelin


// emProBind input file example
//
// format is [ bplist | x3DNAbp | x3DNAparams ]
//
// collection_type = [ EEDR | pulling ]
// pulling_force = {x,y,z}
// base_collection = format:filename
// base_bound_domains = {n1_a:n1_b, n2_a:n2_b, ...} [{}]
// protein_collection = format:filename
// protein_binding_sites = {s1, s2, ...}
// binding_ramp_sampling = N
// base_seqdep_model = model_name [IdealDNA]


#ifndef emDNA_emProBind_h
#define emDNA_emProBind_h


#include <emDNA_Includes.h>
class BpCollection;


// entry point
int main(int argc, char* argv[]);

// collection type
enum class CollectionType : Integer {
    ClampedCollection = 0,
    PullCollection = 1
};

// input parsing function
struct emProBindData {
    CollectionType _coll_type;
    std::string _pulling_force;
    std::string _base_input;
    std::string _base_bound_domains;
    std::string _base_seqdep_model;
    std::string _protein_input;
    std::string _protein_binding_sites;
    Size _ramp_sampling;
};
emProBindData parse_input(int argc, char* argv[]);

// bp collection creation function
BpCollection create_bp_collection(const std::string& input_string);

// inspector message function
void print_inspector_output(const BpCollection& base_collection,
                            const BpCollection& protein_collection,
                            const std::vector<Size>& binding_indices);

// protein_binding_sites key parser function
std::vector<Size> parse_protein_binding_sites(const std::string& s);


#endif  // emDNA_emProBind_h
