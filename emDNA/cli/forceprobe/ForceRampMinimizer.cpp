// ForceRampMinimizer class
// Nicolas Clauvelin


#include <emDNA.h>
#include <ForceRampMinimizer.h>


// initial bp collection accessor/modifier
const BpCollection& ForceRampMinimizer::initial_bp_collection() const {
    return m_initial_conf;
};
void ForceRampMinimizer::set_initial_bp_collection(const BpCollection&
                                                   bp_collection) {
    m_initial_conf = bp_collection;
};


// settings modifiers
void ForceRampMinimizer::set_force_ramp(const ForceRamp& force_ramp) {
    m_force_ramp = force_ramp;
};


// minimization method
bool ForceRampMinimizer::force_ramp_minimization() {

    // minimization data
    std::vector<Real> energy;
    std::vector<Real> gradient_norms_angular;
    std::vector<Real> gradient_norms_distance;
    std::vector<BpCollection> collections;
    std::vector<std::string> return_codes;
    std::vector<Size> iterations;

    // starting configuration
    BpCollection current_conf = m_initial_conf;

    // output files
    setup_output_files();

    // force increment
    const Real F_ini = m_force_ramp._initial_force;
    const Real F_fin = m_force_ramp._final_force;
    const Real delta_f = (F_fin-F_ini)/m_force_ramp._n_increments;
    const Vector3 F_dir = m_force_ramp._force_direction;

    // minimization loop
    // we add 1 to the loop to get the desired starting and final values
    for (Size i=0; i<m_force_ramp._n_increments+1; ++i) {

        // current pulling force
        Real force_magnitude = F_ini + i*delta_f;
        Vector3 pulling_force = force_magnitude*F_dir;

        // interface
        std::shared_ptr<PullingBpCollection> bp_coll_intf_ptr =
        std::make_shared<PullingBpCollection>();
        bp_coll_intf_ptr->set_bp_collection(current_conf);
        bp_coll_intf_ptr->set_pulling_force(pulling_force);

        // progress indicator
        print_progress_indicator(i);

        // callback options
        CallbackOptions callback_opts;
        callback_opts._callback_output = false;

        // minimization
        MinimizationResults minim_res =
        MinimizerAgent::alglib_minimization(bp_coll_intf_ptr,
                                            callback_opts,
                                            nullptr);

        // check return code
        if (minim_res._return_code == "BADCONDS" ||
            minim_res._return_code == "BADGRAD" ||
            minim_res._return_code == "UNKNOWN" ||
            minim_res._return_code == "MAXIT")
            return false;

        // new optimized collection
        current_conf = minim_res._optimized_bp_collection;

        // results - stats
        std::vector<std::string> data;
        data.push_back(EnhancedString::convert_to_string(pulling_force.norm()));
        data.push_back(EnhancedString::convert_to_string(minim_res._final_E));
        data.push_back(EnhancedString::
                       convert_to_string(minim_res._free_dofs_gradient.norm()));
        data.push_back(EnhancedString::
                       convert_to_string(minim_res._n_iterations));
        data.push_back(minim_res._return_code);
        std::stringstream output;
        Parser_x3DNA::fixed_width_vector_formatting(output, data, 8, 15);
        m_pulling_stats_file.open_in_append_mode();
        m_pulling_stats_file.write_stream(output);
        m_pulling_stats_file.close();

        // results - collection
        m_collection_file.open_in_append_mode();
        m_collection_file.
        write_line(format_bp_collection(minim_res._optimized_bp_collection));
        m_collection_file.close();

    };

    // final collection
    OutputFileHandler output_file("emDNA_fp_final.txt",
                                  ".");
    Size n_bps = current_conf.n_of_base_pairs();
    output_file.open();
    for (Size i=0; i<n_bps; ++i) {
        std::stringstream bp_str;
        bp_str << current_conf.base_pair(i);
        bp_str << "\n";
        output_file.write_stream(bp_str);
    };
    output_file.write_line("");
    output_file.close();
    
    // clear the line
    std::cout << "\n";
    std::cout.flush();
    
    return true;

};


// output files setup
void ForceRampMinimizer::setup_output_files() {

    // file names
    m_pulling_stats_file.set_file_name("emDNA_fp_stats.txt");
    m_pulling_stats_file.set_file_path(".");
    m_collection_file.set_file_name("emDNA_fp_confs.txt");
    m_collection_file.set_file_path(".");

    // create empty files
    m_pulling_stats_file.open();
    m_pulling_stats_file.close();
    m_collection_file.open();
    m_collection_file.close();

    // binding stats header
    std::vector<std::string> header;
    header.push_back("FORCE");
    header.push_back("ENERGY");
    header.push_back("|GRAD_ANGLE|");
    header.push_back("|GRAD_DIST|");
    header.push_back("ITERATIONS");
    header.push_back("RETURN_CODE");
    std::stringstream output;
    Parser_x3DNA::fixed_width_vector_formatting(output, header, 8, 15);
    m_pulling_stats_file.open();
    m_pulling_stats_file.write_stream(output);
    m_pulling_stats_file.close();

};


// collection formating
std::string ForceRampMinimizer::
format_bp_collection(const BpCollection& bp_collection) const {

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
void ForceRampMinimizer::print_progress_indicator(Size p) const {
    std::cout << "[";
    for (Size i=0; i<m_force_ramp._n_increments; ++i) {
        if (i <= p)
            std::cout << "|";
        else
            std::cout << " ";
    };
    std::cout << "]\r";
    std::cout.flush();
};
