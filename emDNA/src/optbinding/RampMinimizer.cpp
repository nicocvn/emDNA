// RampMinimizer class
// Nicolas Clauvelin


#include <ProteinBindingRamp.h>
#include <RampMinimizer.h>


// constants
#define CONFS_FILE "emDNA_pb_confs.txt"
#define STATS_FILE "emDNA_pb_stats.txt"
#define FINAL_FILE "emDNA_pb_final.txt"


// initial bp collection accessor/modifier
const BpCollection& RampMinimizer::initial_bp_collection() const {
    return m_initial_conf;
};
void RampMinimizer::set_initial_bp_collection(const BpCollection& bp_collection)
{
    m_initial_conf = bp_collection;
};


// protein bp collection accessor/modifier
const BpCollection& RampMinimizer::protein_bp_collection() const {
    return m_protein_conf;
};
void RampMinimizer::set_protein_bp_collection(const BpCollection&
                                              protein_bp_collection) {
    m_protein_conf = protein_bp_collection;
};


// settings modifiers
void RampMinimizer::set_protein_binding_indices(const std::vector<Size>&
                                                binding_indices) {
    m_binding_indices = binding_indices;
};
void RampMinimizer::set_binding_ramp_sampling(Size binding_ramp_sampling) {
    m_binding_ramp_sampling = binding_ramp_sampling;
};


// output files setup
void RampMinimizer::setup_output_files() {

    // file names
    m_binding_stats_file.set_file_name(STATS_FILE);
    m_binding_stats_file.set_file_path(".");
    m_collection_file.set_file_name(CONFS_FILE);
    m_collection_file.set_file_path(".");

    // create empty files
    m_binding_stats_file.open();
    m_binding_stats_file.close();
    m_collection_file.open();
    m_collection_file.close();

    // binding stats header
    std::stringstream output;
    Parser_x3DNA::fixed_width_vector_formatting(output,
                                                headers_list(),
                                                8, 15);
    m_binding_stats_file.open();
    m_binding_stats_file.write_stream(output);
    m_binding_stats_file.close();

};
std::vector<std::string> RampMinimizer::headers_list() const {
    return {
        "ENERGY", "|GRAD_NORM|", "ITERATIONS", "RETURN_CODE"
    };
};
void RampMinimizer::log_stats(const std::vector<std::string>& stats) {
    std::stringstream output;
    Parser_x3DNA::fixed_width_vector_formatting(output, stats, 8, 15);
    m_binding_stats_file.open_in_append_mode();
    m_binding_stats_file.write_stream(output);
    m_binding_stats_file.close();
};


// collection formating
std::string
RampMinimizer::format_bp_collection(const BpCollection& bp_collection) const {

    // opening
    std::stringstream output;
    output << "{";

    // loop
    Size n_bps = bp_collection.n_of_base_pairs();
    for (Size i=0; i<n_bps-1; ++i)
        output << bp_collection.base_pair(i) << ", ";
    output << bp_collection.base_pair(n_bps-1) << "}";

    return output.str();

};


// progress indicator
void RampMinimizer::print_progress_indicator(Size p) const {
    std::cout << "  " << p+1 << " / " << m_binding_ramp_sampling << "\r";
    std::cout.flush();
};


#undef CONFS_FILE
#undef STATS_FILE
#undef FINAL_FILE

