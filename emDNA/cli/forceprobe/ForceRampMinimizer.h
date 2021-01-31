// ForceRampMinimizer class
// Nicolas Clauvelin


#ifndef emDNA_ForceRampMinimizer_h
#define emDNA_ForceRampMinimizer_h


#include <BpCollection.h>


struct ForceRamp {
    Vector3 _force_direction;
    Real _initial_force;
    Real _final_force;
    Size _n_increments;
};


class ForceRampMinimizer {


public:

    // constructors
    ForceRampMinimizer() = default;
    ForceRampMinimizer(const ForceRampMinimizer& f_ramp_min) = default;
    ~ForceRampMinimizer() = default;

    // copy operator
    ForceRampMinimizer&
    operator=(const ForceRampMinimizer& f_ramp_min) = default;

    // initial bp collection accessor/modifier
    const BpCollection& initial_bp_collection() const;
    void set_initial_bp_collection(const BpCollection& bp_collection);

    // settings modifiers
    void set_force_ramp(const ForceRamp& force_ramp);

    // minimization method
    bool force_ramp_minimization();


private:

    // output files setup
    void setup_output_files();

    // collection formating
    std::string format_bp_collection(const BpCollection& bp_collection) const;

    // progress indicator
    void print_progress_indicator(Size p) const;

    BpCollection m_initial_conf;
    ForceRamp m_force_ramp;

    // output files
    OutputFileHandler m_pulling_stats_file;
    OutputFileHandler m_collection_file;
    
    
};



#endif  // emDNA_ForceRampMinimizer_h
