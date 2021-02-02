// PDBChain class
// Nicolas Clauvelin


#ifndef DNASim_PDBChain_h
#define DNASim_PDBChain_h


#include "pdb/PDBResidue.h"


namespace DNASim {


    class PDBChain {


    public:

        // constructors
        PDBChain() = default;
        PDBChain(const PDBChain& pdb_chain) = default;
        PDBChain(PDBChain&& pdb_chain) = default;
        ~PDBChain() = default;

        // copy and move operators
        PDBChain& operator=(const PDBChain& pdb_chain) = default;
        PDBChain& operator=(PDBChain&& pdb_chain) = default;

        // chain residues management methods
        void add_residue(const PDBResidue& residue);
        void clear();

        // chain identifier accessor/modifier
        const std::string& identifier() const;
        void set_identifier(const std::string& chain_id);

        // chain residue accessor/modifier
        const PDBResidue& chain_residue(Size index) const;
        PDBResidue& chain_residue(Size index);

        // chain residue accessor/modifier by sequence index
        const PDBResidue& chain_residue_from_sequence_index(Size index) const;
        PDBResidue& chain_residue_from_sequence_index(Size index);

        // residue iterators
        std::vector<PDBResidue>::const_iterator residue_begin() const;
        std::vector<PDBResidue>::const_iterator residue_end() const;

        // chain properties accessors
        Size n_atoms() const;
        Size n_residues() const;
        std::string chain_sequence() const;
        bool is_DNA() const;
        

    private:

        std::string m_chain_id;
        std::vector<PDBResidue> m_chain_residues;

    };


}


#endif  // DNASim_PDBChain_h
