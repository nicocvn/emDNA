// ffield_packager tool
// Nicolas Clauvelin


#include <TCLAP_Includes.h>
#include <emDNA.h>
#include <ffield_packager.h>


#define SEQ_DIM 4
#define STEP_DIM 6


// entry point
int main(int argc, char* argv[]) {

    try {

        // command line parsing
        CLData cl_data = parse_command_line(argc, argv);

        // file parsing
        TupleTVector<VectorN> steps_data =
        parse_text_file(cl_data._steps_filename);
        TupleTVector<VectorN> fmat_data =
        parse_text_file(cl_data._fmat_filename);

        // db building
        StepParametersDB steps_db = build_steps_db(cl_data._model_name,
                                                   steps_data);
        ForceConstantsDB fmat_db = build_fmat_db(cl_data._model_name,
                                                 fmat_data);

        // archive
        archive_force_field<CerealBinaryOutput>(cl_data._model_name+".ff",
                                                steps_db, fmat_db);

        // archive
        archive_force_field<CerealJSONOutput>(cl_data._model_name+".json",
                                              steps_db, fmat_db);

        // message
        std::cout << "\n";
        std::cout << "force field packaged in:\n";
        std::cout << "  " + cl_data._model_name + ".ff (binary version)\n";
        std::cout << "  " + cl_data._model_name + ".json (text version)\n";
        std::cout << "\n";

    }

    catch (DNASim_ExitException& e) {
        exit(e._exit_code);
    };

};


// command line parsing option
CLData parse_command_line(int argc, char* argv[]) {

    // command line parser
    TCLAP::CmdLine cmd_line("emDNA_force_field - Nicolas Clauvelin, "
                            "Rutgers University",
                            '=',
                            "");

    // input arguments
    TCLAP::ValueArg<std::string> steps_input("", "intrinsic-steps-input",
                                             "Intrinsic step parameters input "
                                             "file.",
                                             true,
                                             "",
                                             "string");
    TCLAP::ValueArg<std::string> fmat_input("", "force-constants-input",
                                            "Fore constants input file.",
                                            true,
                                            "",
                                            "string");
    TCLAP::ValueArg<std::string> model_name("", "model-name",
                                            "Name for the sequence dependence "
                                            "model",
                                            true,
                                            "",
                                            "string");

    // parsing
    cmd_line.add(steps_input);
    cmd_line.add(fmat_input);
    cmd_line.add(model_name);
    cmd_line.parse(argc, argv);

    // input settings
    CLData cl_data;
    cl_data._steps_filename = steps_input.getValue();
    cl_data._fmat_filename = fmat_input.getValue();
    cl_data._model_name = model_name.getValue();

    return cl_data;

};


//Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
// text file parsing function for data
TupleTVector<VectorN> parse_text_file(const std::string& filename) {

    // file slurping
    std::vector<std::string> file_lines =
    InputFileHandler::read_file(filename, ".");

    // container
    TupleTVector<VectorN> data;

    data.reserve((SEQ_DIM+1)*SEQ_DIM*SEQ_DIM*(SEQ_DIM+1));

    std::set<std::string> tetramers;//Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)


    if(file_lines.size()==16){
        for (Size i=0; i<file_lines.size(); ++i) {

            // string splitting
            std::vector<std::string> tokens;
            bool tk_flag = EnhancedString::tokenize_string(file_lines[i],
                                                        tokens,
                                                        '=');
            DS_ASSERT(tk_flag,
                    "wrong format for sequence dependence file");

            // data
            const char *seq = tokens[0].c_str();
            StepSequence step_seq(Base::base_symbol_from_char(seq[0]),
                                Base::base_symbol_from_char(seq[1]));
            VectorN flatten_vals(tokens[1]);

            // tuple
            const BaseSymbol n1 = step_seq.first_base();
            const BaseSymbol n2 = step_seq.last_base();

            char bases[SEQ_DIM+1] = {'A', 'C', 'T', 'G', 'x'};
            for(Size i=0; i<SEQ_DIM+1; ++i){
                for(Size j=0; j<SEQ_DIM+1; ++j){
                    data.push_back(TupleT<VectorN>(Base::base_symbol_from_char(bases[i]), n1, n2,
                                                    Base::base_symbol_from_char(bases[j]), flatten_vals));
                };
            };
        };

    }
    //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    else if (file_lines.size()==400){
         for (Size i=0; i<file_lines.size(); ++i) {

            // string splitting
            std::vector<std::string> tokens;
            bool tk_flag = EnhancedString::tokenize_string(file_lines[i],
                                                        tokens,
                                                        '=');
            DS_ASSERT(tk_flag,
                    "wrong format for sequence dependence file");

            // data
            const char *seq = tokens[0].c_str();
            tetramers.insert(seq);//Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
            TetramerSequence tetra_seq(Base::base_symbol_from_char(seq[0]),
                                Base::base_symbol_from_char(seq[1]),
                                Base::base_symbol_from_char(seq[2]),
                                Base::base_symbol_from_char(seq[3]));
            VectorN flatten_vals(tokens[1]);

            // tuple
            const BaseSymbol n1 = tetra_seq.first_base();
            const BaseSymbol n2 = tetra_seq.second_base();
            const BaseSymbol n3 = tetra_seq.third_base();
            const BaseSymbol n4 = tetra_seq.fourth_base();
            data.push_back(TupleT<VectorN>(n1, n2, n3, n4, flatten_vals));
        };

        //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
        //If set less than 400, then some tetramer was added more than once
        DS_ASSERT(tetramers.size() == 400,
                "The input file " + filename + " contains duplicate tetramers");

    }
    //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    else {
        DS_ASSERT(false,
              "the input file " + filename + " is missing data; "
              "the file has to contain the 16 (or 400) sequence combinations");
    }

    return data;

};

// step parameters db building function
//Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
StepParametersDB build_steps_db(const std::string& model_name,
                                const TupleTVector<VectorN>& data) {

    // data re-organization
    StepParameters step_data[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1];

    // loop
    for (Size i=0; i<data.size(); ++i) {

        // size checking
        DS_ASSERT(std::get<4>(data[i]).size() == STEP_DIM,
                  "wrong dimension for step paramteters data");

        // step parameters
        step_data[(Size)std::get<0>(data[i])][(Size)std::get<1>(data[i])]
                    [(Size)std::get<2>(data[i])][(Size)std::get<3>(data[i])] =
        StepParameters(std::get<4>(data[i]));

    };

    return StepParametersDB(model_name, step_data);

};

// force constants db building function
//Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
ForceConstantsDB build_fmat_db(const std::string& model_name,
                               const TupleTVector<VectorN>& data) {

    // data re-organization
    MatrixN fmat_data[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1];

    // loop
    for (Size i=0; i<data.size(); ++i) {

        // size checking
        DS_ASSERT(std::get<4>(data[i]).size() == STEP_DIM*STEP_DIM,
                  "wrong dimension for force constant matrix data");

        // matrix
        MatrixN fmat(STEP_DIM, STEP_DIM);
        fmat << std::get<4>(data[i]).pack_values();
        
        // step parameters
        fmat_data[(Size)std::get<0>(data[i])][(Size)std::get<1>(data[i])] 
                    [(Size)std::get<2>(data[i])][(Size)std::get<3>(data[i])] =
        fmat;
        
    };
    
    return ForceConstantsDB(model_name, fmat_data);
    
};


#undef SEQ_DIM
#undef STEP_DIM
