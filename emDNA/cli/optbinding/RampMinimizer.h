// RampMinimizer class
// Nicolas Clauvelin


// virtual implementation of the ramp minimizer; this class should be derived
// and adapted to a specific end conditions case
// the only method to adapt is binding_ramp_minimization()


#ifndef emDNA_RampMinimizer_h
#define emDNA_RampMinimizer_h


#include <emDNA.h>


class RampMinimizer {


public:

    // constructors
    RampMinimizer() = default;
    RampMinimizer(const RampMinimizer& ramp_min) = default;
    RampMinimizer(RampMinimizer&& ramp_min) = default;
    virtual ~RampMinimizer() = default;

    // copy and move operators
    RampMinimizer& operator=(const RampMinimizer& ramp_min) = default;
    RampMinimizer& operator=(RampMinimizer&& ramp_min) = default;

    // initial bp collection accessor/modifier
    const BpCollection& initial_bp_collection() const;
    void set_initial_bp_collection(const BpCollection& bp_collection);

    // protein bp collection accessor/modifier
    const BpCollection& protein_bp_collection() const;
    void set_protein_bp_collection(const BpCollection& protein_bp_collection);

    // settings modifiers
    void set_protein_binding_indices(const std::vector<Size>& binding_indices);
    void set_binding_ramp_sampling(Size binding_ramp_sampling);

    // minimization method
    virtual bool binding_ramp_minimization() =0;


protected:

    // output files method
    void setup_output_files();
    virtual std::vector<std::string> headers_list() const;
    virtual void log_stats(const std::vector<std::string>& stats);

    // collection formating
    std::string format_bp_collection(const BpCollection& bp_collection) const;

    // progress indicator
    void print_progress_indicator(Size p) const;

    BpCollection m_initial_conf;            // base configuration
    BpCollection m_protein_conf;            // protein configuration
    std::vector<Size> m_binding_indices;    // index of the binding step
    Size m_binding_ramp_sampling;           // binding ramp sampling

    // output files
    OutputFileHandler m_binding_stats_file;
    OutputFileHandler m_collection_file;


};


#endif  // emDNA_RampMinimizer_h
