// SequenceDependenceModels header
// Nicolas Clauvelin


#include "dna/SequenceDepenceModels.h"


namespace DNASim {


// sequence-dependence model structure
SequenceDependenceModelData::
SequenceDependenceModelData(const std::string step_params[16],
                            const std::string force_constants[16],
                            const std::string& description) :
    _step_parameters_data(),
    _force_constants_data(),
    _description(description) {
        for (Size i=0; i<16; ++i) {
            _step_parameters_data[i] = step_params[i];
            _force_constants_data[i] = force_constants[i];
        };
};


}
