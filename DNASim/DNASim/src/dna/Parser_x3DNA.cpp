// Parser_x3DNA class implementation file
// Nicolas Clauvelin


#include <EnhancedString.h>
#include <InputFileHandler.h>
#include <Vector3.h>
#include <StepParameters.h>
#include <Triad.h>
#include <Sequence.h>
#include <Parser_x3DNA.h>


namespace DNASim {


    // x3DNA step parameters file parsing method
    // format description:
    // X-Y  ..  ..  ..  ..  ..  ..  rho1  rho2  rho3  theta1  theta2  theta3
    Data_x3DNA<StepParameters>
    Parser_x3DNA::read_step_parameters_file(const std::string& filename) {

        // file reading
        std::vector<std::string> file_lines(InputFileHandler::
                                            read_file(filename, "."));

        // data container
        std::vector<std::string> bp_sequences;
        std::vector<StepParameters> step_parameters;

        // loop over all lines
        for (Size i=0, end=file_lines.size(); i<end; ++i) {

            // discard line if starting with a # (comment) or of zero length
            // lines containing a # are also discarded
            if (file_lines[i][0] == '#' ||          // starts with a #
                file_lines[i].size() == 0 ||        // zero length
                std::find(file_lines[i].begin(),    // contains a #
                          file_lines[i].end(), '#') != file_lines[i].end())
                continue;

            // line processing
            // split the string using whitespaces
            std::vector<std::string> tokens =
            EnhancedString::tokenize_string(file_lines[i], ' ');

            // first token is the step sequence
            // we need to remove the '-'
            bp_sequences.push_back(EnhancedString::
                                   remove_character(tokens[0], '-'));

            // step parameters are given by tokens 8,..,13
            VectorN prms(6, FLOAT_INIT);
            for (Size j=0; j<6; ++j)
#define x3DNA_PRMS_SHIFT 7
                prms[j] =
                EnhancedString::
                convert_from_string<Real>(tokens[j+x3DNA_PRMS_SHIFT]);
#undef x3DNA_PRMS_SHIFT

            step_parameters.push_back(StepParameters(VectorN( {
                prms[3], prms[4], prms[5],
                prms[0],prms[1],prms[2]
            } )));

        };

        // sequence reconstruction
        // we pick the first nucleotide of each bp sequence
        std::string full_sequence;
        for (Size i=0; i<bp_sequences.size(); ++i)
            full_sequence += std::string(&(bp_sequences[i][0]), 1);
        
        // remove first set of step parameters
        // the first line is with zero parameters to get the correct sequence
        // and should be discarded
        step_parameters.erase(step_parameters.begin());
        
        // data container
        Data_x3DNA<StepParameters> parsed_data;
        parsed_data._data = step_parameters;
        parsed_data._sequence = full_sequence;
        
        return parsed_data;

    };


    // x3DNA base-pair frames file parsing method
    // format description:
    //    ...     3 T-A ...
    //    -0.5714539118  1.3371475710   7.8270802256 # origin
    //    0.7349719905  0.6768517980  -0.0410830457 # x-axis
    //    -0.6342840908  0.6647963617  -0.3946257589 # y-axis
    //    -0.2397912952  0.3160972018   0.9179230326 # z-axis
    Data_x3DNA<Triad>
    Parser_x3DNA::read_base_pairs_file(const std::string& filename) {

        // file reading
        std::vector<std::string> file_lines(InputFileHandler::
                                            read_file(filename, "."));

        // loop over all lines for groups
        std::vector<Size> group_indexes;
        for (Size i=0, end=file_lines.size(); i<end; ++i) {

            // discard line if starting with a # (comment) or of zero length
            if (file_lines[i][0] == '#' || file_lines[i].size() == 0)
                continue;

            // found the first line of a group
            if (file_lines[i][0] == '.')
                group_indexes.push_back(i);

        };

        // data container
        std::string sequence;
        std::vector<Triad> base_pairs;

        // group processing
        for (Size i=0; i<group_indexes.size(); ++i) {

            Size idx = group_indexes[i];

            // sequence letter
            sequence += Parser_x3DNA::extract_sequence_letter(file_lines[idx]);

            // base pair
            Triad bp;
            bp.set_origin(Parser_x3DNA::
                          extract_vector_from_bp_line(file_lines[idx+1]));
            bp.set_axis(I,
                        Parser_x3DNA::
                        extract_vector_from_bp_line(file_lines[idx+2]));
            bp.set_axis(J,
                        Parser_x3DNA::
                        extract_vector_from_bp_line(file_lines[idx+3]));
            bp.set_axis(K,
                        Parser_x3DNA::
                        extract_vector_from_bp_line(file_lines[idx+4]));
            bp.orthogonalize();
            base_pairs.push_back(bp);
            
        };
        
        // data container
        Data_x3DNA<Triad> parsed_data;
        parsed_data._data = base_pairs;
        parsed_data._sequence = sequence;
        
        return parsed_data;

    };


    // x3DNA step parameters formatting method
    void Parser_x3DNA::
    create_step_parameters_output(std::stringstream& output,
                                  const Data_x3DNA<StepParameters>& data) {

#define STEP_PARAMS_PREC 5
#define STEP_PARAMS_FIELD_WIDTH 11

        // sizing header
        const Size n_bps = data._data.size()+1;
        const Size n_params = data._data.size();
        output << n_bps << " # base pairs\n";
        output << n_params << " # ***local base-pair & step parameters***\n";

        // header
        std::vector<std::string> header;
        header.push_back("#");
        header.push_back("Shear");
        header.push_back("Stretch");
        header.push_back("Stagger");
        header.push_back("Buckle");
        header.push_back("Prop-Tw");
        header.push_back("Opening");
        header.push_back("Shift");
        header.push_back("Slide");
        header.push_back("Rise");
        header.push_back("Tilt");
        header.push_back("Roll");
        header.push_back("Twist");
        fixed_width_vector_formatting(output,
                                      header,
                                      STEP_PARAMS_PREC,
                                      STEP_PARAMS_FIELD_WIDTH);

        // build base-pair sequences
        Size n_steps = data._data.size();
        std::vector<std::string> sequences;
        for (Size i=0; i<n_steps; ++i) {
            const char n1 = data._sequence[i];
            const BaseSymbol n2s =
            Base::complementary_base(Base::base_symbol_from_char(n1));
            const char n2 = Base::str(n2s)[0];
            sequences.push_back(std::string(&n1,1)+"-"+std::string(&n2,1));
        };
        const char n1 = data._sequence[n_steps];
        const BaseSymbol n2s =
        Base::complementary_base(Base::base_symbol_from_char(n1));
        const char n2 = Base::str(n2s)[0];
        sequences.push_back(std::string(&n1,1)+"-"+std::string(&n2,1));

        // loop over all steps
        for (Size i=0; i<data._data.size()+1; ++i)
            if (i == 0) {
                x3DNA_format_bp_step_params(output,
                                            StepParameters(),
                                            sequences[i]);
            }
            else {
                x3DNA_format_bp_step_params(output,
                                            data._data[i-1],
                                            sequences[i]);
            };

#undef STEP_PARAMS_PREC
#undef STEP_PARAMS_FIELD_WIDTH

    };


    // x3DNA base pairs formatting method
    void Parser_x3DNA::create_base_pairs_output(std::stringstream& output,
                                                const Data_x3DNA<Triad>& data) {

        // sizing header
        const Size n_bps = data._data.size();
        output << n_bps << " base pairs\n";

        // loop over all base pairs
        for (Size i=0; i<data._data.size(); ++i) {

            // bp sequence
            std::string bp_seq;
            const char n1 = data._sequence[i];
            const BaseSymbol n2s =
            Base::complementary_base(Base::base_symbol_from_char(n1));
            const char n2 = Base::str(n2s)[0];
            bp_seq = std::string(&n1);
            bp_seq += "-";
            bp_seq += std::string(&n2);
            
            // formatting
            x3DNA_format_base_pair(output, data._data[i], i+1, bp_seq);
            
        };

    };


    // private parsing methods
    std::string Parser_x3DNA::extract_sequence_letter(const std::string& s) {

        // string tokenizing with whitespace
        std::vector<std::string> tokens =
        EnhancedString::tokenize_string(s, ' ');

        // take the first character of the third token
#define TK_IDX 2
        return std::string(&tokens[TK_IDX][0], 1);
#undef TK_IDX

    };
    Vector3 Parser_x3DNA::extract_vector_from_bp_line(const std::string& s) {

        // string tokenizing with whitespace
        std::vector<std::string> tokens =
        EnhancedString::tokenize_string(s, ' ');

        // vector
        return Vector3(EnhancedString::convert_from_string<Real>(tokens[0]),
                       EnhancedString::convert_from_string<Real>(tokens[1]),
                       EnhancedString::convert_from_string<Real>(tokens[2]));

    };


    // x3DNA bp step params formatting method
#define STEP_PARAMS_PREC 5
#define STEP_PARAMS_FIELD_WIDTH 11
    // sequence is expected as "X-Y"
    void Parser_x3DNA::x3DNA_format_bp_step_params(std::stringstream& output,
                                                   const StepParameters& p,
                                                   const std::string& seq) {

        // line elements
        std::vector<std::string> str_elements;

        // sequence
        str_elements.push_back(seq);

        // zero parameters (account for base step parameters)
#define x3DNA_N_OF_BASE_STEP_PARAMS 6
        for (Size i=0; i<x3DNA_N_OF_BASE_STEP_PARAMS; ++i)
            str_elements.
            push_back(fixed_width_convert_to_string(FLOAT_INIT,
                                                    STEP_PARAMS_PREC));
#undef x3DNA_N_OF_BASE_STEP_PARAMS

        // bp step params values
        str_elements.
        push_back(fixed_width_convert_to_string(p[SHIFT], STEP_PARAMS_PREC));
        str_elements.
        push_back(fixed_width_convert_to_string(p[SLIDE], STEP_PARAMS_PREC));
        str_elements.
        push_back(fixed_width_convert_to_string(p[RISE], STEP_PARAMS_PREC));
        str_elements.
        push_back(fixed_width_convert_to_string(p[TILT], STEP_PARAMS_PREC));
        str_elements.
        push_back(fixed_width_convert_to_string(p[ROLL], STEP_PARAMS_PREC));
        str_elements.
        push_back(fixed_width_convert_to_string(p[TWIST], STEP_PARAMS_PREC));

        // formatted string
        fixed_width_vector_formatting(output, str_elements,
                                      STEP_PARAMS_PREC,
                                      STEP_PARAMS_FIELD_WIDTH);
        
    };
#undef STEP_PARAMS_PREC
#undef STEP_PARAMS_FIELD_WIDTH


    // x3DNA base pair formatting method
#define BP_FIELD_WIDTH 16
    void Parser_x3DNA::x3DNA_format_base_pair(std::stringstream& output,
                                              const Triad& bp,
                                              const Size& bp_index,
                                              const std::string& bp_seq) {

        // base pair sequence
        output << "... " << bp_index << " " << bp_seq << " ...\n";

        // base pair origin
        fixed_width_vector_formatting(output, bp.origin().std_vector(),
                                      REAL_WIDTH, BP_FIELD_WIDTH);

        // base pair axes
        fixed_width_vector_formatting(output, bp.axis(I).std_vector(),
                                      REAL_WIDTH, BP_FIELD_WIDTH);
        fixed_width_vector_formatting(output, bp.axis(J).std_vector(),
                                      REAL_WIDTH, BP_FIELD_WIDTH);
        fixed_width_vector_formatting(output, bp.axis(K).std_vector(),
                                      REAL_WIDTH, BP_FIELD_WIDTH);
        
    };
#undef BP_FIELD_WIDTH



}

