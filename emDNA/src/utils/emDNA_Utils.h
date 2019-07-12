// emDNA_Utils class
// Nicolas Clauvelin


// collection of useful functions


#ifndef emDNA_emDNA_Utils_h
#define emDNA_emDNA_Utils_h


#include <emDNA_Includes.h>


namespace emDNA_Utils {


    // frozen steps parsing functions
    // first function:
    // parse a list "x1:x2, y1:y2, ..." and returns a sorted vector of ranges
    // [x1:x2], [y1:y2], ...
    // second function:
    // same as first but parse a list of the form "{...}"
    std::vector<SizePair> parse_frozen_steps_list(const std::string& s);
    std::vector<SizePair> parse_frozen_steps_vector(const std::string& s);

    // parse minim settings vector
    // parse a vector of the form {max-iterations,dx,f,g,max-step-size} and
    // returns a AlglibMinSettings structure
    AlglibMinSettings parse_minim_settings(const std::string& s);

}


#endif  // emDNA_emDNA_Utils_h
