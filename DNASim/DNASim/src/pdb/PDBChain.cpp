// PDBChain class
// Nicolas Clauvelin


#include <PDBChain.h>


namespace DNASim {


    // chain residues management methods
    void PDBChain::add_residue(const PDBResidue& residue) {
        return m_chain_residues.push_back(residue);
    };
    void PDBChain::clear() {
        m_chain_residues.clear();
    };


    // chain identifier accessor/modifier
    const std::string& PDBChain::identifier() const {
        return m_chain_id;
    };
    void PDBChain::set_identifier(const std::string& chain_id) {
        m_chain_id = chain_id;
    };


    // chain residue accessor/modifier
    const PDBResidue& PDBChain::chain_residue(Size index) const {
        return m_chain_residues[index];
    };
    PDBResidue& PDBChain::chain_residue(Size index) {
        return m_chain_residues[index];
    };


    // chain residue accessor/modifier by sequence index
    const PDBResidue& PDBChain::
    chain_residue_from_sequence_index(Size index) const {
        auto it = std::find_if(m_chain_residues.begin(), m_chain_residues.end(),
                               [&index](const PDBResidue& res) {
                                   return (index ==
                                           res.residue_sequence_index());
                               });
        return (*it);
    };
    PDBResidue& PDBChain::chain_residue_from_sequence_index(Size index) {
        auto it = std::find_if(m_chain_residues.begin(), m_chain_residues.end(),
                               [&index](const PDBResidue& res) {
                                   return (index ==
                                           res.residue_sequence_index());
                               });
        return (*it);
    };


    // residue iterators
    std::vector<PDBResidue>::const_iterator PDBChain::residue_begin() const {
        return m_chain_residues.begin();
    };
    std::vector<PDBResidue>::const_iterator PDBChain::residue_end() const {
        return m_chain_residues.end();
    };


    // chain properties accessors
    Size PDBChain::n_atoms() const {
        Size natoms = Size(0);
        std::for_each(m_chain_residues.begin(),
                      m_chain_residues.end(),
                      [&natoms](const PDBResidue& res){
                          natoms += res.n_atoms();
                      });
        return natoms;
    };
    Size PDBChain::n_residues() const {
        return m_chain_residues.size();
    };
    std::string PDBChain::chain_sequence() const {
        std::stringstream chain_seq;
        for (const PDBResidue& res : m_chain_residues)
            chain_seq << res.name() << ":";
        std::string seq(chain_seq.str());
        return seq.substr(0, seq.size()-1);
    };
    bool PDBChain::is_DNA() const {
        std::string first_res = m_chain_residues.front().name();
        if (first_res == "A" || first_res == "DA")
            return true;
        else if (first_res == "C" || first_res == "DC")
            return true;
        else if (first_res == "T" || first_res == "DT")
            return true;
        else if (first_res == "G" || first_res == "DG")
            return true;
        else
            return false;
    };


}
