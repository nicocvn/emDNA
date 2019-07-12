// RampMinimizerEEDR class
// Nicolas Clauvelin


#include <Minimizer_Alglib.h>
#include <MinimizerAgent.h>
#include <ClampedBpCollection.h>
#include <ProteinBindingRamp.h>
#include <RampMinimizerEEDR.h>


// constants
#define CONFS_FILE "emDNA_pb_confs.txt"
#define STATS_FILE "emDNA_pb_stats.txt"
#define FINAL_FILE "emDNA_pb_final.txt"


// minimization method
bool RampMinimizerEEDR::binding_ramp_minimization() {

    // number of proteins
    Size n_proteins = m_binding_indices.size();

    // binding domain size and corresponding frozen domains
    const Size binding_domain_size = m_protein_conf.n_of_bp_steps();
    std::vector<SizePair> frozen_domains =
    m_initial_conf.frozen_steps_domains();
    for (Size i=0; i<n_proteins; ++i) {
        frozen_domains.push_back(SizePair(m_binding_indices[i],
                                          m_binding_indices[i]
                                          +binding_domain_size-1));
    };

    // protein step parameters
    const std::vector<BpStepParams> protein_prms(m_protein_conf.
                                                 bp_step_params().begin(),
                                                 m_protein_conf.
                                                 bp_step_params().end());

    // base configuration step parameters for the binding domains
    std::vector<std::vector<BpStepParams>> base_prms;
    for (Size k=0; k<n_proteins; ++k) {
        std::vector<BpStepParams> prms;
        for (Size i=0; i<binding_domain_size; ++i)
            prms.
            push_back(m_initial_conf.bp_step_params(m_binding_indices[k]+i));
        base_prms.push_back(prms);
    };

    // protein binding ramps
    std::vector<ProteinBindingRamp> binding_ramps(m_binding_indices.size());
    for (Size k=0; k<n_proteins; ++k) {
        binding_ramps[k].set_base_configuration_parameters(base_prms[k]);
        binding_ramps[k].set_protein_configuration_parameters(protein_prms);
    };

    // minimization data
    std::vector<Real> energy;
    std::vector<Real> gradient_norms_angular;
    std::vector<Real> gradient_norms_distance;
    std::vector<BpCollection> collections;
    std::vector<std::string> return_codes;
    std::vector<Size> iterations;

    // starting configuration
    BpCollection current_conf = m_initial_conf;
    current_conf.set_frozen_steps_domains(frozen_domains);

    // original first base pair
    const BasePair first_bp = m_initial_conf.base_pair(0);

    // output files
    setup_output_files();

    // minimization loop
    Real delta = 1./(Real(m_binding_ramp_sampling));
    for (Size i=0; i<m_binding_ramp_sampling; ++i) {

        // ramp step parameters
        std::vector<std::vector<BpStepParams>> ramps_prms;
        for (Size k=0; k<n_proteins; ++k) {
            ramps_prms.push_back(binding_ramps[k].
                                 generate_step_parameters(delta*(i+1)));
        };

        // copy all step parameters from the current solution
        std::vector<BpStepParams>
        guess_prms(current_conf.bp_step_params().begin(),
                   current_conf.bp_step_params().end());

        // copy the binding steps
        for (Size k=0; k<n_proteins; ++k) {
            std::copy(ramps_prms[k].begin(), ramps_prms[k].end(),
                      guess_prms.begin()+m_binding_indices[k]);
        };

        // before updating the end base pair we look for the last non-frozen
        // step
        // we look for the last free step:
        const SizePair& last_frozen_size_pair =
        current_conf.frozen_steps_domains().back();
        const Size& last_frozen_step_index = last_frozen_size_pair.second;
        Size free_step_index = 0;
        if (last_frozen_step_index < current_conf.n_of_bp_steps()-1)
            // the last step is free
            free_step_index = current_conf.n_of_bp_steps()-1;
        else
            // we take the step right before the frozen domain
            free_step_index = last_frozen_size_pair.first-1;

        // collection update
        std::vector<BasePair> base_pairs =
        BpStepParams::rebuild_bps(guess_prms, first_bp);
        for (Size i=free_step_index+1; i<m_initial_conf.n_of_base_pairs(); ++i)
            base_pairs[i] = m_initial_conf.base_pair(i);
        base_pairs[free_step_index+1] =
        m_initial_conf.base_pair(free_step_index+1);
        current_conf.update_collection_from_base_pairs(base_pairs);

        // interface
        std::shared_ptr<BpCollection_Interface> bp_coll_intf_ptr =
        std::make_shared<ClampedBpCollection>();
        bp_coll_intf_ptr->set_bp_collection(current_conf);

        // callback options
        CallbackOptions callback_opts;
        callback_opts._callback_output = false;

        // progress indicator
        print_progress_indicator(i);

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
        current_conf.
        update_collection_from_base_pairs(minim_res.
                                          _optimized_bp_collection.
                                          base_pairs());

        // results - stats
        std::vector<std::string> data;
        data.push_back(EnhancedString::convert_to_string(minim_res._final_E));
        data.push_back(EnhancedString::
                       convert_to_string(minim_res._free_dofs_gradient.norm()));
        data.push_back(EnhancedString::
                       convert_to_string(minim_res._n_iterations));
        data.push_back(minim_res._return_code);
        log_stats(data);

        // results - collection
        m_collection_file.open_in_append_mode();
        m_collection_file.
        write_line(format_bp_collection(minim_res._optimized_bp_collection));
        m_collection_file.close();

    };

    // final collection
    OutputFileHandler output_file(FINAL_FILE, ".");
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


// constants
#undef CONFS_FILE
#undef STATS_FILE
#undef FINAL_FILE
