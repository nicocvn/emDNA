// PDBStructure class
// Nicolas Clauvelin


#ifndef DNASim_PDBStructure_h
#define DNASim_PDBStructure_h


#include "pdb/PDBChain.h"


namespace DNASim {


    class PDBStructure {


    public:

        // constructors
        PDBStructure() = default;
        PDBStructure(const PDBStructure& pdb_structure) = default;
        PDBStructure(PDBStructure&& pdb_structure) = default;
        ~PDBStructure() = default;

        // copy and move operators
        PDBStructure& operator=(const PDBStructure& pdb_structure) = default;
        PDBStructure& operator=(PDBStructure&& pdb_structure) = default;

        // structure chain management methods
        void add_chain(const PDBChain& chain);
        void clear();

        // structure name accessor/modifier
        const std::string& name() const;
        void set_name(const std::string& structure_name);

        // structure chain accessor/modifier
        const PDBChain& structure_chain(Size index) const;
        PDBChain& structure_chain(Size index);

        // structure chain accessor/modifier by chain identifier
        const PDBChain& structure_chain(const std::string& chain_id) const;
        PDBChain& structure_chain(const std::string& chain_id);

        // chain iterators
        std::vector<PDBChain>::const_iterator chain_begin() const;
        std::vector<PDBChain>::const_iterator chain_end() const;

        // chain properties accessors
        Size n_atoms() const;
        Size n_residues() const;
        Size n_chains() const;


    private:

        std::string m_structure_name;
        std::vector<PDBChain> m_structure_chains;
        std::map<std::string, Size> m_chain_id_table;


    };


}


#endif  // DNASim_PDBStructure_h
