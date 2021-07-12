// ForceConstantsDB class
// Nicolas Clauvelin


#include "file_io/EnhancedString.h"
#include "maths/VectorN.h"
#include "dna/Sequence.h"
#include "dna/ForceConstantsDB.h"

// sequence dependence models
#include "dna/SequenceDepenceModels.h"

// tetramer dependence models
#include "dna/TetramerDepenceModels.h"


namespace DNASim {


    // macro for sequence dimension
    #define SEQ_DIM 4


    // class constructor with initialization from a model name
    ForceConstantsDB::ForceConstantsDB(const std::string& model_name) :
    m_force_constants(),
    m_model_name() {
        init_with_model(model_name);
    };


    // class constructor with full initialization (model name + db data)
    ForceConstantsDB::ForceConstantsDB(const std::string& model_name,
                                       const MatrixN (*db_data)[SEQ_DIM]) :
    m_force_constants(),
    m_model_name(model_name) {

        // size checking
        for (Size i=0; i<SEQ_DIM; ++i) {
            const Size sz_chk = (Size)(sizeof(db_data[i])/sizeof(MatrixN));
            DS_ASSERT(sz_chk==SEQ_DIM,
                      "wrong size for db_data");
        };

        // initialization
        for (Size i=0; i<(SEQ_DIM+1); ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                for (Size k=0; k<SEQ_DIM; ++k) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
                    for (Size l=0; l<(SEQ_DIM+1); ++l) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
                        m_force_constants[i][j][k][l] = db_data[j][k]; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)

    };


    // class constructor with full initialization (model name + db data)
    ForceConstantsDB::ForceConstantsDB(const std::string& model_name,
                                       const MatrixN (db_data)[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1]) :
    m_force_constants(),
    m_model_name(model_name) {

        // size checking
        for (Size i=0; i<SEQ_DIM+1; ++i) {
            const Size sz_chk = (Size)(sizeof(db_data[i])/sizeof(MatrixN));
            DS_ASSERT(sz_chk==SEQ_DIM*SEQ_DIM*(SEQ_DIM+1),
                      "wrong size for db_data");
        };

        // initialization
        for (Size i=0; i<(SEQ_DIM+1); ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                for (Size k=0; k<SEQ_DIM; ++k) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
                    for (Size l=0; l<(SEQ_DIM+1); ++l) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
                        m_force_constants[i][j][k][l] = db_data[i][j][k][l]; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)

    };

    // DB init methods
    void ForceConstantsDB::init_with_model(const std::string& model_name) {

       
        // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
        auto it1 = SequenceDependenceModelList.find(model_name);
        auto it2 = TetramerDependenceModelList.find(model_name);
        

        if (it1 != SequenceDependenceModelList.end()){
            init_data_from_list_of_string(it1->second._force_constants_data);
        }
        else if (it2 != TetramerDependenceModelList.end()) {
            init_tetrameric_data_from_list_of_string(it2->second._force_constants_data);
        }
        else {
            DS_ASSERT(false,
            "non-existing sequence-dependence model name: " + model_name);
        }
        
        // model name
        m_model_name = model_name;

    };


    // sequence dependent force constants
    const MatrixN&
    ForceConstantsDB::force_constants(const StepSequence& step_seq) const {
        const Size i = static_cast<Size>(step_seq.first_base());
        const Size j = static_cast<Size>(step_seq.last_base());
        return m_force_constants[4][i][j][4];
    };

    // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    // tetrameric dependent force constants
    const MatrixN&
    ForceConstantsDB::force_constants(const TetramerSequence& tetra_seq) const {
        const Size i = static_cast<Size>(tetra_seq.first_base());
        const Size j = static_cast<Size>(tetra_seq.second_base());
        const Size k = static_cast<Size>(tetra_seq.third_base());
        const Size l = static_cast<Size>(tetra_seq.fourth_base());
        return m_force_constants[i][j][k][l];
    };


    // model name accessor
    const std::string& ForceConstantsDB::model_name() const {
        return m_model_name;
    };


    // IdealDNA init method
    void ForceConstantsDB::
    init_data_from_list_of_string(const std::string data[SEQ_DIM*SEQ_DIM]) {

        // filling two dimensional array
        for (Size i=0; i<SEQ_DIM*SEQ_DIM; ++i) {

            // string splitting
            std::vector<std::string> tokens;
            bool chk = EnhancedString::tokenize_string(data[i], tokens, '=');

            // checking
            DS_ASSERT(chk && tokens[0].size()==2,
                      "wrong format in intrinsic step parameters data\n"
                      "string = " + data[i]);

            // step sequence
            const char *seq = tokens[0].c_str();
            StepSequence step_seq(Base::base_symbol_from_char(seq[0]),
                                  Base::base_symbol_from_char(seq[1]));

            // data
            VectorN flatten_mat(tokens[1]);
#define STEP_DIM 6
            MatrixN fmat(STEP_DIM);
#undef STEP_DIM
            fmat << flatten_mat.pack_values();

            // set up
            const Size idx1 = static_cast<Size>(step_seq.first_base());
            const Size idx2 = static_cast<Size>(step_seq.last_base());

            for (Size j=0; j<(SEQ_DIM+1); ++j) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
                for (Size k=0; k<(SEQ_DIM+1); ++k) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
                    m_force_constants[j][idx1][idx2][k] = fmat; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
            
        };
        
    };

    //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    void ForceConstantsDB::
    init_tetrameric_data_from_list_of_string(const std::string data[(SEQ_DIM+1)*SEQ_DIM*SEQ_DIM*(SEQ_DIM+1)]) {

        // filling two dimensional array
        for (Size i=0; i<(SEQ_DIM+1)*SEQ_DIM*SEQ_DIM*(SEQ_DIM+1); ++i) {

            // string splitting
            std::vector<std::string> tokens;
            bool chk = EnhancedString::tokenize_string(data[i], tokens, '=');

            // checking
            DS_ASSERT(chk && tokens[0].size()==4,
                      "wrong format in intrinsic step parameters data\n"
                      "string = " + data[i]);

            // step sequence
            const char *seq = tokens[0].c_str();
            TetramerSequence tetra_seq(Base::base_symbol_from_char(seq[0]),
                                  Base::base_symbol_from_char(seq[1]),
                                  Base::base_symbol_from_char(seq[2]),
                                  Base::base_symbol_from_char(seq[3]));

            // data
            VectorN flatten_mat(tokens[1]);
#define STEP_DIM 6
            MatrixN fmat(STEP_DIM);
#undef STEP_DIM
            fmat << flatten_mat.pack_values();

            // set up
            const Size idx1 = static_cast<Size>(tetra_seq.first_base());
            const Size idx2 = static_cast<Size>(tetra_seq.second_base());
            const Size idx3 = static_cast<Size>(tetra_seq.third_base());
            const Size idx4 = static_cast<Size>(tetra_seq.fourth_base());
            
            m_force_constants[idx1][idx2][idx3][idx4] = fmat; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
            
        };
        
    };



    #undef SEQ_DIM


}

