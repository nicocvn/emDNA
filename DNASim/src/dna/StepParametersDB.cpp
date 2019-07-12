// StepParametersDB class
// Nicolas Clauvelin


#include <EnhancedString.h>
#include <Sequence.h>
#include <StepParametersDB.h>

// sequence dependence models
#include <SequenceDepenceModels.h>


namespace DNASim {


    // macro for sequence dimension
    #define SEQ_DIM 4


    // class constructor with initialization from a model name
    StepParametersDB::StepParametersDB(const std::string& model_name) :
    m_step_parameters(),
    m_model_name() {
        init_with_model(model_name);
    };


    // class constructor with full initialization (model name and db data)
    // this constructor is designed for serialization purpose only
    StepParametersDB::StepParametersDB(const std::string& model_name,
                                       const StepParameters
                                       (*db_data)[SEQ_DIM]) :
    m_step_parameters(),
    m_model_name(model_name) {

        // size checking
        for (Size i=0; i<SEQ_DIM; ++i) {
            const Size sz_chk =
            (Size)(sizeof(db_data[i])/sizeof(StepParameters));
            DS_ASSERT(sz_chk==SEQ_DIM,
                      "wrong size for db_data");
        };

        // initialization
        for (Size i=0; i<SEQ_DIM; ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                m_step_parameters[i][j] = db_data[i][j];

    };


    // DB init methods
    void StepParametersDB::init_with_model(const std::string& model_name) {

        // look for model name
        auto it = SequenceDependenceModelList.find(model_name);
        DS_ASSERT(it != SequenceDependenceModelList.end(),
                  "non-existing sequence-dependence model name: " + model_name);

        // initialization
        init_data_from_list_of_string(it->second._step_parameters_data);

        // model name
        m_model_name = model_name;

    };


    // sequence dependent step parameters
    const StepParameters&
    StepParametersDB::intrinsic_bp_step_params(const StepSequence& step_seq)
    const {
        const Size i = static_cast<Size>(step_seq.first_base());
        const Size j = static_cast<Size>(step_seq.last_base());
        return m_step_parameters[i][j];
    };


    // model name accessor/modifier
    const std::string& StepParametersDB::model_name() const {
        return m_model_name;
    };
    void StepParametersDB::set_model_name(const std::string& name) {
        m_model_name = name;
    };


    // IdealDNA init method
    void StepParametersDB::
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
            StepParameters p(tokens[1]);

            // set up
            const Size idx1 = static_cast<Size>(step_seq.first_base());
            const Size idx2 = static_cast<Size>(step_seq.last_base());
            m_step_parameters[idx1][idx2] = p;

        };

    };


    #undef SEQ_DIM


}
