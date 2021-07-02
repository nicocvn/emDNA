// ForceConstantsDB class
// Nicolas Clauvelin


// class for sequence-dependent force constants management


#ifndef MC_DNA_ForceConstantsDB_h
#define MC_DNA_ForceConstantsDB_h


#include "maths/MatrixN.h"


namespace DNASim {


    class StepSequence;
    class TetramerSequence; // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)


    // macro for sequence dimension
    #define SEQ_DIM 4


    class ForceConstantsDB {


    public:

        // constructors
        ForceConstantsDB() = default;
        ForceConstantsDB(const ForceConstantsDB& force_constants_db)  = default;
        ForceConstantsDB(ForceConstantsDB&& force_constants_db)  = default;
        ForceConstantsDB(const std::string& model_name);
        ForceConstantsDB(const std::string& model_name,
                         const MatrixN (*db_data)[SEQ_DIM]);
        ForceConstantsDB(const std::string& model_name,
                         const MatrixN (db_data)[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1]);
        ~ForceConstantsDB()  = default;

        // copy and move operators
        ForceConstantsDB& operator=(const ForceConstantsDB&
                                    force_constants_db) = default;
        ForceConstantsDB& operator=(ForceConstantsDB&&
                                    force_constants_db) = default;

        // DB init methods
        void init_with_model(const std::string& model_name);

        // sequence dependent step parameters
        const MatrixN& force_constants(const StepSequence& step_seq) const;
        
        //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
        // tetramer dependent force constants
        const MatrixN& force_constants(const TetramerSequence& tetra_seq) const;


        // model name accessor
        const std::string& model_name() const;


    private:
        // model init methods
        void init_data_from_list_of_string(const std::string
                                           data[SEQ_DIM*SEQ_DIM]);

        //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
        // model init methods
        void init_tetrameric_data_from_list_of_string(const std::string
                                           data[(SEQ_DIM+1)*SEQ_DIM*SEQ_DIM*(SEQ_DIM+1)]);

        // sequence-dependent data
        MatrixN m_force_constants[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1]; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)

        // model name
        std::string m_model_name;
        
        
    };


    #undef SEQ_DIM


}


#endif  // MC_DNA_ForceConstantsDB_h
