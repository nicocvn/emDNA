// StepParametersDB class
// Nicolas Clauvelin


// class for sequence-dependent intrinsic step parameters management


#ifndef MC_DNA_StepParametersDB_h
#define MC_DNA_StepParametersDB_h


#include <StepParameters.h>


namespace DNASim {


    class StepSequence;


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

        // model name accessor/modifier
        const std::string& model_name() const;
        void set_model_name(const std::string& name);


    private:

        // model init methods
        void init_data_from_list_of_string(const std::string
                                           data[SEQ_DIM*SEQ_DIM]);

        // sequence-dependent data
        StepParameters m_step_parameters[SEQ_DIM][SEQ_DIM];

        // model name
        std::string m_model_name;
        
        
    };


    #undef SEQ_DIM


}


#endif  // MC_DNA_StepParametersDB_h
