// TetramerDepenceModels header
// Zoe Wefers (McGill University, June 2021, DIMACS REU)


// list of available depence-dependence models
//
// the list of models is implemented as a std::map<model_name, model_data>


#ifndef emDNA_TetramerDepenceModels_h
#define emDNA_TetramerDepenceModels_h


#include "DNASim_Includes.h"


namespace DNASim {


    // list of tetramer dependence model - step parameters
    #include "dna/StepParameters_Cohen2017_tet.h"

    // list of tetramer dependence model - force constants
    #include "dna/ForceConstants_Cohen2017_tet.h"


    // tetramer-dependence model structure
    struct TetramerDependenceModelData {
        std::string _step_parameters_data[400];
        std::string _force_constants_data[400];
        std::string _description;
        TetramerDependenceModelData(const std::string step_params[400],
                                    const std::string force_constants[400],
                                    const std::string& description);
    };


    // list of implemented models
    static std::map<std::string,TetramerDependenceModelData>
    TetramerDependenceModelList {

        // TetramericModel
        { "Cohen2017_tet",
            TetramerDependenceModelData(StepParameters_Cohen2017_tet,
                                        ForceConstants_Cohen2017_tet,
                                        "Tetrameric force field from Cohen (2017)")
        }

    };


}


#endif  // emDNA_TetramerDepenceModels_h
