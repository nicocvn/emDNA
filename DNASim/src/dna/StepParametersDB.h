// StepParametersDB class
// Nicolas Clauvelin


// class for sequence-dependent intrinsic step parameters management


#ifndef MC_DNA_StepParametersDB_h
#define MC_DNA_StepParametersDB_h


#include "dna/StepParameters.h"


namespace DNASim {


    class StepSequence;
    class TetramerSequence; //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)


    // macro for sequence dimension
    #define SEQ_DIM 4


    class StepParametersDB {


    public:

        // constructors
        // the constructor with *db_data is designed for serialization purpose
        // only
        StepParametersDB() = default;
        StepParametersDB(const StepParametersDB& step_parameters_db) = default;
        StepParametersDB(StepParametersDB&& step_parameters_db) = default;
        StepParametersDB(const std::string& model_name);
        StepParametersDB(const std::string& model_name,
                         const StepParameters (*db_data)[SEQ_DIM]);
        StepParametersDB(const std::string& model_name,
                         const StepParameters (db_data)[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1]);
        ~StepParametersDB() = default;

        // copy and move operators
        StepParametersDB& operator=(const StepParametersDB&
                                    step_parameters_db) = default;
        StepParametersDB& operator=(StepParametersDB&&
                                    step_parameters_db) = default;

        // DB init methods
        void init_with_model(const std::string& model_name);

        // sequence dependent step parameters
        const StepParameters&
        intrinsic_bp_step_params(const StepSequence& step_seq) const;

        // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
        const StepParameters&
        intrinsic_bp_step_params(const TetramerSequence& step_seq) const;

        // model name accessor/modifier
        const std::string& model_name() const;
        void set_model_name(const std::string& name);


    private:

        // model init methods
        void init_data_from_list_of_string(const std::string
                                           data[SEQ_DIM*SEQ_DIM]);

        // model init methods
        void init_tetrameric_data_from_list_of_string(const std::string
                                           data[(SEQ_DIM+1)*SEQ_DIM*SEQ_DIM*(SEQ_DIM+1)]); //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)

        // sequence-dependent data
        StepParameters m_step_parameters[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1]; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)

        // model name
        std::string m_model_name;
        
        
    };


    #undef SEQ_DIM


}


#endif  // MC_DNA_StepParametersDB_h
