// Sequence classes
// Nicolas Clauvelin


// multiple classes for the management of DNA sequence


#ifndef MC_DNA_Sequence_h
#define MC_DNA_Sequence_h


#include "DNASim_Includes.h"


namespace DNASim {


    // nucleotide type enum
    enum class BaseSymbol : Size {
        A = 0, C = 1, G = 2, T = 3, x = 4 // "x" added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    };


    // Base class
    class Base {


    public:

        // static methods for manipulating BaseSymbol
        static std::string str(const BaseSymbol& base);
        static BaseSymbol base_symbol_from_char(const char s);
        static BaseSymbol complementary_base(const BaseSymbol& base);


    private:

        // private constructors and operators
        Base() = delete;
        Base(const Base& n) = delete;
        Base(Base&& n) = delete;
        Base& operator=(const Base& n) = delete;
        Base& operator=(Base&& n) = delete;
        ~Base() = delete;


    };


    // StepSequence class
    class StepSequence {


    public:

        // constructors
        StepSequence() = default;
        StepSequence(const BaseSymbol& n1, const BaseSymbol& n2);
        StepSequence(const StepSequence& step_sequence) = default;
        StepSequence(StepSequence&& step_sequence) = default;
        ~StepSequence() = default;

        // copy and move operators
        StepSequence& operator=(const StepSequence& step_sequence) = default;
        StepSequence& operator=(StepSequence&& step_sequence) = default;

        // base accessors
        const BaseSymbol& first_base() const;
        const BaseSymbol& last_base() const;


    private:

        std::pair<BaseSymbol,BaseSymbol> m_bases;


    };

    //Added by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    //TetramerSequence class
    class TetramerSequence{
    
    public:
        //constructors
        TetramerSequence() = default;
        TetramerSequence(const BaseSymbol& n1, const BaseSymbol& n2, const BaseSymbol& n3, const BaseSymbol& n4);
        TetramerSequence(const TetramerSequence& tetra_sequence) = default;
        TetramerSequence(TetramerSequence&& tetra_sequence) = default;
        ~TetramerSequence() = default;

        // copy and move operators (not sure what these are for)
        TetramerSequence& operator=(const TetramerSequence& tetra_sequence) = default;
        TetramerSequence& operator=(TetramerSequence&& tetra_sequence) = default;

        // base accessors
        const BaseSymbol& first_base() const;
        const BaseSymbol& second_base() const;
        const BaseSymbol& third_base() const;
        const BaseSymbol& fourth_base() const;

    private:
        std::tuple<BaseSymbol, BaseSymbol, BaseSymbol, BaseSymbol> m_bases;

    };




    // Sequence class
    class Sequence {


    public:

        // constructors
        Sequence() = default;
        Sequence(const std::string& linear_sequence);
        Sequence(const Sequence& sequence) = default;
        Sequence(Sequence&& sequence) = default;
        ~Sequence() = default;

        // copy and move operators
        Sequence& operator=(const Sequence& sequence) = default;
        Sequence& operator=(Sequence&& sequence) = default;

        // step sequence accessor
        const StepSequence& step_sequence(const Size& step_index) const;

        // linear sequence rebuilding method
        std::string linear_sequence() const;


    private:

        std::vector<StepSequence> m_step_sequences;


    };


}


#endif  // MC_DNA_Sequence_h
