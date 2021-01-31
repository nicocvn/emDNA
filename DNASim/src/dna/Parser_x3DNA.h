// Parser_x3DNA class header file
// Nicolas Clauvelin


#ifndef DNASim_Parser_x3DNA_h
#define DNASim_Parser_x3DNA_h


#include "DNASim_Includes.h"


namespace DNASim {


    class Vector3;
    class StepParameters;
    class Triad;


    // data container with sequence
    template <class DataType>
    struct Data_x3DNA {
        std::vector<DataType> _data;
        std::string _sequence;
    };


    class Parser_x3DNA {


    public:

        // static x3DNA reading methods
        static Data_x3DNA<StepParameters>
        read_step_parameters_file(const std::string& filename);
        static Data_x3DNA<Triad>
        read_base_pairs_file(const std::string& filename);

        // static x3DNA writing methods
        static void
        create_step_parameters_output(std::stringstream& output,
                                      const Data_x3DNA<StepParameters>& data);
        static void
        create_base_pairs_output(std::stringstream& output,
                                 const Data_x3DNA<Triad>& data);


        // templatized formatting methods
        template <class T>
        static void fixed_width_vector_formatting(std::stringstream& output,
                                                  const std::vector<T>& vec,
                                                  Size precision_width,
                                                  Size fixed_width) {
            output << std::fixed << std::setprecision((Integer)precision_width);
            for (Size i=0; i<vec.size()-1; ++i)
                output << std::left << std::setw((Integer)fixed_width)
                    << vec[i];
            output << vec[vec.size()-1] << "\n";
        };
        template <class T>
        static std::string fixed_width_convert_to_string(const T& x,
                                                         Size precision_width =
                                                         REAL_WIDTH) {
            std::stringstream ss;
            ss << std::fixed << std::setprecision((Integer)precision_width)
                << x;
            return std::string(ss.str());
        };


    private:

        // private parsing methods
        static std::string extract_sequence_letter(const std::string& s);
        static Vector3 extract_vector_from_bp_line(const std::string& s);

        // private formatting methods
        static void x3DNA_format_bp_step_params(std::stringstream& output,
                                                const StepParameters& p,
                                                const std::string& seq);
        static void x3DNA_format_base_pair(std::stringstream& output,
                                           const Triad& bp,
                                           const Size& bp_index,
                                           const std::string& bp_seq);

        // constructors and operators
        Parser_x3DNA() = delete;
        Parser_x3DNA(const Parser_x3DNA& parser) = delete;
        Parser_x3DNA(Parser_x3DNA&& parser) = delete;
        Parser_x3DNA& operator=(const Parser_x3DNA& parser) = delete;
        Parser_x3DNA& operator=(Parser_x3DNA&& parser) = delete;


    };


}


#endif  // DNASim_Parser_x3DNA_h
