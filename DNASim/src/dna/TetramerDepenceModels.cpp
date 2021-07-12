// TetramerDepenceModels struct
// Zoe Wefers (McGill University, June 2021, DIMACS REU)


#include "dna/TetramerDepenceModels.h"


namespace DNASim {


// Tetramer-dependence model structure
TetramerDependenceModelData::
TetramerDependenceModelData(const std::string step_params[400],
                            const std::string force_constants[400],
                            const std::string& description) :
    _step_parameters_data(),
    _force_constants_data(),
    _description(description) {
        for (Size i=0; i<400; ++i) {
            _step_parameters_data[i] = step_params[i];
            _force_constants_data[i] = force_constants[i];
        };
};


}