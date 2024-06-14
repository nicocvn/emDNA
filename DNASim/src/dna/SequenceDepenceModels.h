// SequenceDependenceModels header
// Nicolas Clauvelin


// list of available sequence-dependence models
//
// the list of models is implemented as a std::map<model_name, model_data>


#ifndef emDNA_SequenceDepenceModels_h
#define emDNA_SequenceDepenceModels_h


#include "DNASim_Includes.h"


namespace DNASim {


    // list of sequence dependence model - step parameters
    #include "dna/StepParameters_IdealDNA.h"
    #include "dna/StepParameters_IdealDNA_304.h"
    #include "dna/StepParameters_AnisoDNA.h"
    #include "dna/StepParameters_AnisoDNA_304.h"
    #include "dna/StepParameters_Olson1998.h"
    #include "dna/StepParameters_Olson1998Legacy.h"
    #include "dna/StepParameters_Cohen2017_dim.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)

    // list of sequence dependence model - force constants
    #include "dna/ForceConstants_IdealDNA.h"
    #include "dna/ForceConstants_AnisoDNA.h"
    #include "dna/ForceConstants_Olson1998.h"
    #include "dna/ForceConstants_Cohen2017_dim.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)


    // sequence-dependence model structure
    struct SequenceDependenceModelData {
        std::string _step_parameters_data[16];
        std::string _force_constants_data[16];
        std::string _description;
        SequenceDependenceModelData(const std::string step_params[16],
                                    const std::string force_constants[16],
                                    const std::string& description);
    };


    // list of implemented models
    static std::map<std::string,SequenceDependenceModelData>
    SequenceDependenceModelList {

        // IdealDNA
        { "IdealDNA",
            SequenceDependenceModelData(StepParameters_IdealDNA,
                                        ForceConstants_IdealDNA,
                                        "ideal B-DNA like isotropic force "
                                        "field")
        },

        // IdealDNA
        { "IdealDNA_304",
            SequenceDependenceModelData(StepParameters_IdealDNA_304,
                                        ForceConstants_IdealDNA,
                                        "ideal B-DNA like isotropic force "
                                        "field (twist=34.3)")
        },

        // AnisoDNA
        { "AnisoDNA",
            SequenceDependenceModelData(StepParameters_AnisoDNA,
                                        ForceConstants_AnisoDNA,
                                        "ideal B-DNA like anisotropic force "
                                        "field")
        },

        // AnisoDNA_304
        { "AnisoDNA_304",
            SequenceDependenceModelData(StepParameters_AnisoDNA_304,
                                        ForceConstants_AnisoDNA,
                                        "ideal B-DNA like anisotropic force "
                                        "field (twist=34.3)")
        },

        // Olson1998
        { "Olson1998",
            SequenceDependenceModelData(StepParameters_Olson1998,
                                        ForceConstants_Olson1998,
                                        "force field from W. K. Olson 98 paper")
        },

        // Olson1998 legacy
        { "Olson1998Legacy",
            SequenceDependenceModelData(StepParameters_Olson1998Legacy,
                                        ForceConstants_Olson1998,
                                        "force field from W. K. Olson 98 "
                                        "paper (legacy)")
        },

        // Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
        // Cohen2017_dim
        { "Cohen2017_dim",
            SequenceDependenceModelData(StepParameters_Cohen2017_dim,
                                        ForceConstants_Cohen2017_dim,
                                        "Dimeric force field from Cohen (2017)")
        },

    };


}


#endif  // emDNA_SequenceDepenceModels_h
