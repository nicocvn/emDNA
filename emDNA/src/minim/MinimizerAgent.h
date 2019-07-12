// MinimizerAgent
// Nicolas Clauvelin


#ifndef emDNA_MinimizerAgent_h
#define emDNA_MinimizerAgent_h


#include <emDNA_Includes.h>
#include <BpCollection.h>
#include <Minimizer_Alglib.h>
#include <BpCollection_Interface.h>


// minimization results data structure
struct MinimizationResults {

    // optimized bp colleciton
    BpCollection _optimized_bp_collection;

    // initial and final energy values
    Real _initial_E;
    Real _final_E;

    // list of contributions to the energy (final)
    std::vector<EnergyContrib> _energy_contribs;

    // elastic energy gradient
    VectorN _elastic_gradient;

    // free dofs gradient (gradient used in the mininmization)
    VectorN _free_dofs_gradient;

    // minimization outcome parameters
    Size _n_iterations;
    std::string _return_code;
    std::string _minim_desc;

};


// typedef
using AlgligMinSettings_Ptr = std::unique_ptr<const AlglibMinSettings>;


class MinimizerAgent {


public:

    // alglib minimization method
    static MinimizationResults
    alglib_minimization(const std::shared_ptr<BpCollection_Interface>&
                        bp_intf_ptr,
                        const CallbackOptions& callback_opts,
                        const AlgligMinSettings_Ptr& settings);


private:

    // default minim settings
    static AlglibMinSettings default_minim_settings();

    // private constructors and copy operator
    MinimizerAgent() = delete;
    MinimizerAgent(const MinimizerAgent& minim_agent) = delete;
    MinimizerAgent(MinimizerAgent&& minim_agent) = delete;
    MinimizerAgent& operator=(const MinimizerAgent& minim_agent) = delete;
    MinimizerAgent& operator=(MinimizerAgent&& minim_agent) = delete;
    ~MinimizerAgent() = delete;


};


#endif  // emDNA_MinimizerAgent_h
