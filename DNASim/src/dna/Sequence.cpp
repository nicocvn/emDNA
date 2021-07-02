// Sequence classes
// Nicolas Clauvelin


#include "dna/Sequence.h"


namespace DNASim {


    /*** Base class ***/

    // static methods
    std::string Base::str(const BaseSymbol& base) {
        // the following array has to be ordered as the BaseSymbol enum
        const std::string str[5] = {"A", "C", "G", "T", "x"};
        return str[static_cast<Size>(base)];
    };
    BaseSymbol Base::base_symbol_from_char(const char s) {
        if (s == 'A')
            return BaseSymbol::A;
        else if (s == 'C')
            return BaseSymbol::C;
        else if (s == 'G')
            return BaseSymbol::G;
        else if (s == 'T')
            return BaseSymbol::T;
        else if (s == 'x')  // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
            return BaseSymbol::x;
        else
            DS_ASSERT(false, "wrong base symbol:\n"+std::string(&s));
    };
    BaseSymbol Base::complementary_base(const BaseSymbol& base) {
        if (base == BaseSymbol::A)
            return BaseSymbol::T;
        else if (base == BaseSymbol::T)
            return BaseSymbol::A;
        else if (base == BaseSymbol::C)
            return BaseSymbol::G;
        else if (base == BaseSymbol::G)
            return BaseSymbol::C;
        else if (base == BaseSymbol::x) //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
            return BaseSymbol::x;
        else
            DS_ASSERT(false, "wrong base symbol:\n"+str(base));
    };


    /*** StepSequence class ***/

    // constructor with full initialization
    StepSequence::StepSequence(const BaseSymbol& base_1,
                               const BaseSymbol& base_2) :
    m_bases(base_1, base_2) {};


    // base accessors
    const BaseSymbol& StepSequence::first_base() const {
        return m_bases.first;
    };
    const BaseSymbol& StepSequence::last_base() const {
        return m_bases.second;
    };



    // Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    /*** TetramerSequence class ***/
    TetramerSequence::TetramerSequence(const BaseSymbol& base_1, 
                                        const BaseSymbol& base_2,
                                        const BaseSymbol& base_3, 
                                        const BaseSymbol& base_4) :
    m_bases(base_1, base_2, base_3, base_4) {};

    //base acessors
    const BaseSymbol& TetramerSequence::first_base() const {
        return std::get<0>(m_bases);
    };
    const BaseSymbol& TetramerSequence::second_base() const {
        return std::get<1>(m_bases);
    };
    const BaseSymbol& TetramerSequence::third_base() const {
        return std::get<2>(m_bases);
    };
    const BaseSymbol& TetramerSequence::fourth_base() const {
        return std::get<3>(m_bases);
    };

    /*** Sequence class ***/

    // constructor with initialization
    // size is reduced by one because the number of steps is the number of bases
    // minus one
    Sequence::Sequence(const std::string& linear_sequence) {

        // sequence length
        const Size& n_bases = linear_sequence.size();

        // char array
        const char *sequence_chars = linear_sequence.c_str();

        if (linear_sequence.size() != 0) {

            // building step sequences
            m_step_sequences = std::vector<StepSequence>(n_bases-1,
                                                         StepSequence());
            for (Size i=0; i<n_bases-1; ++i) {
                m_step_sequences[i] =
                StepSequence(Base::base_symbol_from_char(sequence_chars[i]),
                             Base::base_symbol_from_char(sequence_chars[i+1]));
            };

        };

    };


    // step sequence accessor
    const StepSequence& Sequence::step_sequence(const Size& step_index) const {
        return m_step_sequences[step_index];
    };


    // linear sequence
    std::string Sequence::linear_sequence() const {

        if (m_step_sequences.size() == 0)
            return std::string();

        // string container
        std::string lin_sequence;

        // loop over all step sequences
        for (const StepSequence& seq : m_step_sequences) {
            lin_sequence += Base::str(seq.first_base());
        };

        // add the terminal base to get the full sequence
        lin_sequence += Base::str(m_step_sequences.back().last_base());

        return lin_sequence;

    };


}

