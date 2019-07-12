// emDNA_Utils class
// Nicolas Clauvelin


#include <emDNA_Utils.h>


namespace emDNA_Utils {

    // frozen steps parsing functions
    std::vector<SizePair> parse_frozen_steps_list(const std::string& s) {

        // check if string is empty (i.e., no frozen steps)
        if (s.size() == 0)
            return std::vector<SizePair>();

        // container
        std::vector<SizePair> frozen_parts;

        // tokenizing for ranges
        std::vector<std::string> range_tokens =
        EnhancedString::tokenize_string(s, ',');

        // loop over the ranges
        for (Size i=0; i<range_tokens.size(); ++i) {

            // tokenizing
            std::vector<std::string> tokens =
            EnhancedString::tokenize_string(range_tokens[i], ':');

            // check format
            DS_ASSERT(tokens.size()==2, "badly formatted frozen steps domain");

            // convert
            // -1 because of the numbering convention
            Size n1 = (Size)EnhancedString::
            convert_from_string<Integer>(tokens[0]);
            Size n2 = (Size)EnhancedString::
            convert_from_string<Integer>(tokens[1]);
            DS_ASSERT(n1 != Size(0) && n2 != Size(0),
                      "steps numbering starts at 1 for frozen steps domain");
            n1 -= 1;
            n2 -= 1;

            // check ordering
            DS_ASSERT(n1<=n2,
                      "wrong range specification for frozen steps domain");

            // inclusive range
            SizePair frozen_steps(static_cast<Size>(n1),
                                  static_cast<Size>(n2));

            frozen_parts.push_back(frozen_steps);

        };

        // sorting
        std::sort(frozen_parts.begin(), frozen_parts.end(),
                  [](const SizePair& sz_pair1, const SizePair& sz_pair2){
                      return sz_pair1.first < sz_pair2.first;
                  });
        
        return frozen_parts;

    };
    std::vector<SizePair> parse_frozen_steps_vector(const std::string& s) {

        // string cleaning
        EnhancedString clean_str(s);
        clean_str.erase_character('{');
        clean_str.erase_character('}');

        // parsing
        return parse_frozen_steps_list(std::string(clean_str));

    };


    // parse minim settings vector
    AlglibMinSettings parse_minim_settings(const std::string& s) {

        // vector parsing
        VectorN settings(s);
        DS_ASSERT(settings.size()==5,
                  "wrong specifications for the minim settings; five "
                  "values are expected");

        // structure
        AlglibMinSettings minim_settings;
        minim_settings._max_iterations = static_cast<Size>(settings[0]);
        minim_settings._threshold_dx = settings[1];
        minim_settings._threshold_f = settings[2];
        minim_settings._threshold_g = settings[3];
        minim_settings._max_step_size = settings[4];
        
        return minim_settings;

    };

}
