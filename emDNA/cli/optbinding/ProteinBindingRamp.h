// ProteinBindingRamp class
// Nicolas Clauvelin


#ifndef emDNA_ProteinBindingRamp_h
#define emDNA_ProteinBindingRamp_h


#include <BpStepParams.h>


class ProteinBindingRamp {


public:

    // constructors
    ProteinBindingRamp();
    ProteinBindingRamp(const ProteinBindingRamp& protein_binding_ramp);
    ~ProteinBindingRamp();

    // copy operator
    ProteinBindingRamp& operator=(const ProteinBindingRamp&
                                  protein_binding_ramp);

    // setting methods
    void set_base_configuration_parameters(const std::vector<BpStepParams>&
                                           base_conf_parameters);
    void set_protein_configuration_parameters(const std::vector<BpStepParams>&
                                              protein_conf_parameters);

    // binding ramp method
    std::vector<BpStepParams> generate_step_parameters(Real ramp_value) const;


private:

    std::vector<BpStepParams> m_base_conf;
    std::vector<BpStepParams> m_protein_conf;


};


#endif  // emDNA_ProteinBindingRamp_h
