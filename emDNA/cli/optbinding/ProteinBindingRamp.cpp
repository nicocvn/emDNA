// ProteinBindingRamp class
// Nicolas Clauvelin


#include <ProteinBindingRamp.h>


// class default constructor
ProteinBindingRamp::ProteinBindingRamp() :
m_base_conf(),
m_protein_conf() {};


// class constructor by copy
ProteinBindingRamp::ProteinBindingRamp(const ProteinBindingRamp&
                                       protein_binding_ramp) :
m_base_conf(protein_binding_ramp.m_base_conf),
m_protein_conf(protein_binding_ramp.m_protein_conf) {};


// class destructor
ProteinBindingRamp::~ProteinBindingRamp() {};


// copy operator
ProteinBindingRamp& ProteinBindingRamp::operator=(const ProteinBindingRamp&
                                                  protein_binding_ramp) {
    m_base_conf = protein_binding_ramp.m_base_conf;
    m_protein_conf = protein_binding_ramp.m_protein_conf;
    return *this;
};


// setting methods
void ProteinBindingRamp::
set_base_configuration_parameters(const std::vector<BpStepParams>&
                                  base_conf_parameters) {
    m_base_conf = base_conf_parameters;
};
void ProteinBindingRamp::
set_protein_configuration_parameters(const std::vector<BpStepParams>&
                                     protein_conf_parameters) {
    m_protein_conf = protein_conf_parameters;
};


// binding ramp method
std::vector<BpStepParams>
ProteinBindingRamp::generate_step_parameters(Real ramp_value) const {

    // set of step parameters
    std::vector<BpStepParams> params(m_base_conf.size(), BpStepParams());

    // coefficients
    const Real c1 = Real(1)-ramp_value;
    const Real c2 = ramp_value;

    // ramping
    for (Size i=0; i<m_base_conf.size(); ++i)
        params[i] =
        BpStepParams(c1*m_base_conf[i].inline_vector() +
                     c2*m_protein_conf[i].inline_vector());

    return params;

};

