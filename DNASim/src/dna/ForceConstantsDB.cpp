// ForceConstantsDB class
// Nicolas Clauvelin


#include <EnhancedString.h>
#include <VectorN.h>
#include <Sequence.h>
#include <ForceConstantsDB.h>

// sequence dependence models
#include <SequenceDepenceModels.h>


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
        for (Size i=0; i<SEQ_DIM; ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                m_force_constants[i][j] = db_data[i][j];

    };


    // DB init methods
    void ForceConstantsDB::init_with_model(const std::string& model_name) {

        // look for model name
        auto it = SequenceDependenceModelList.find(model_name);
        DS_ASSERT(it != SequenceDependenceModelList.end(),
                  "non-existing sequence-dependence model name: " + model_name);

        // initialization
        init_data_from_list_of_string(it->second._force_constants_data);

        // model name
        m_model_name = model_name;

    };


    // sequence dependent step parameters
    const MatrixN&
    ForceConstantsDB::force_constants(const StepSequence& step_seq) const {
        const Size i = static_cast<Size>(step_seq.first_base());
        const Size j = static_cast<Size>(step_seq.last_base());
        return m_force_constants[i][j];
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
            m_force_constants[idx1][idx2] = fmat;
            
        };
        
    };


    #undef SEQ_DIM


}

