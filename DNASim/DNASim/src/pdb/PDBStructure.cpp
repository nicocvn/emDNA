// PDBStructure class
// Nicolas Clauvelin


#include <PDBStructure.h>


namespace DNASim {


    // structure chain management methods
    void PDBStructure::add_chain(const PDBChain& chain) {
        m_structure_chains.push_back(chain);
        m_chain_id_table.
        insert(std::pair<std::string,Size>(chain.identifier(),
                                           m_structure_chains.size()-1));
    };
    void PDBStructure::clear() {
        m_structure_chains.clear();
        m_chain_id_table.clear();
    };


    // structure chain accessor/modifier by chain identifier
    const PDBChain& PDBStructure::
    structure_chain(const std::string& chain_id) const {
        auto it = m_chain_id_table.find(chain_id);
        return m_structure_chains[it->second];
    };
    PDBChain& PDBStructure::structure_chain(const std::string& chain_id) {
        auto it = m_chain_id_table.find(chain_id);
        return m_structure_chains[it->second];
    };


    // structure name accessor/modifier
    const std::string& PDBStructure::name() const {
        return m_structure_name;
    };
    void PDBStructure::set_name(const std::string& structure_name) {
        m_structure_name = structure_name;
    };


    // structure chain accessor/modifier
    const PDBChain& PDBStructure::structure_chain(Size index) const {
        return m_structure_chains[index];
    };
    PDBChain& PDBStructure::structure_chain(Size index) {
        return m_structure_chains[index];
    };


    // chain iterators
    std::vector<PDBChain>::const_iterator PDBStructure::chain_begin() const {
        return m_structure_chains.begin();
    };
    std::vector<PDBChain>::const_iterator PDBStructure::chain_end() const {
        return m_structure_chains.end();
    };


    // chain properties accessors
    Size PDBStructure::n_atoms() const {
        Size natoms = Size(0);
        std::for_each(m_structure_chains.begin(),
                      m_structure_chains.end(),
                      [&natoms](const PDBChain& chain){
                          natoms += chain.n_atoms();
                      });
        return natoms;
    };
    Size PDBStructure::n_residues() const {
        Size nres = Size(0);
        std::for_each(m_structure_chains.begin(),
                      m_structure_chains.end(),
                      [&nres](const PDBChain& chain){
                          nres += chain.n_residues();
                      });
        return nres;
    };
    Size PDBStructure::n_chains() const {
        return m_structure_chains.size();
    };


}
