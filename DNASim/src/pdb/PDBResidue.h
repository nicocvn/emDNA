// PDBResidue class
// Nicolas Clauvelin


#ifndef DNASim_PDBResidue_h
#define DNASim_PDBResidue_h


#include "pdb/PDBAtom.h"


namespace DNASim  {


    class PDBAtomRecord;


    class PDBResidue {


    public:

        // constructors
        PDBResidue() = default;
        PDBResidue(const PDBResidue& pdb_residue) = default;
        PDBResidue(PDBResidue&& pdb_residue) = default;
        ~PDBResidue() = default;

        // copy and move operators
        PDBResidue& operator=(const PDBResidue& pdb_residue) = default;
        PDBResidue& operator=(PDBResidue&& pdb_residue) = default;

        // residue atom management methods
        void add_atom(const PDBAtom& atom);
        void clear();

        // residue name accessor/modifier
        const std::string& name() const;
        void set_name(const std::string& residue_name);

        // residue serial accessor/modifier
        const Integer& residue_sequence_index() const;
        void set_residue_sequence_index(const Integer& index);

        // residue atom accessor/modifier
        const PDBAtom& residue_atom(Size index) const;
        PDBAtom& residue_atom(Size index);

        // atom finding method
        const PDBAtom& find_atom(const std::string& atom_name) const;

        // residue properties accessors
        Size n_atoms() const;

        // instantiation from a PDBAtomRecord
        static PDBResidue create_from_record(const PDBAtomRecord& pdb_record);


    private:

        Integer m_residue_serial;
        std::string m_residue_name;
        std::vector<PDBAtom> m_residue_atoms;


    };


}


#endif  // DNASim_PDBResidue_h

//Alanine                     ALA                         A
//Arginine                    ARG                         R
//Asparagine                  ASN                         N
//Aspartic acid               ASP                         D
//ASP/ASN ambiguous           ASX                         B
//Cysteine                    CYS                         C
//Glutamine                   GLN                         Q
//Glutamic acid               GLU                         E
//GLU/GLN ambiguous           GLX                         Z
//Glycine                     GLY                         G
//Histidine                   HIS                         H
//Isoleucine                  ILE                         I
//Leucine                     LEU                         L
//Lysine                      LYS                         K
//Methionine                  MET                         M
//Phenylalanine               PHE                         F
//Proline                     PRO                         P
//Serine                      SER                         S
//Threonine                   THR                         T
//Tryptophan                  TRP                         W
//Tyrosine                    TYR                         Y
//Unknown                     UNK
//Valine                      VAL                         V
//
