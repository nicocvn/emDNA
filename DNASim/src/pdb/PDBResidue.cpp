// PDBResidue class
// Nicolas Clauvelin


#include "pdb/PDBAtomRecord.h"
#include "pdb/PDBResidue.h"


namespace DNASim {


    // residue serial accessor/modifier
    const Integer& PDBResidue::residue_sequence_index() const {
        return m_residue_serial;
    };
    void PDBResidue::set_residue_sequence_index(const Integer& index) {
        m_residue_serial = index;
    };


    // residue atom management methods
    void PDBResidue::add_atom(const PDBAtom& atom) {
        m_residue_atoms.push_back(atom);
    };
    void PDBResidue::clear() {
        m_residue_atoms.clear();
    };


    // residue name accessor/modifier
    const std::string& PDBResidue::name() const {
        return m_residue_name;
    };
    void PDBResidue::set_name(const std::string& residue_name) {
        m_residue_name = residue_name;
    };


    // residue atom accessor/modifier
    const PDBAtom& PDBResidue::residue_atom(Size index) const {
        return m_residue_atoms[index];
    };
    PDBAtom& PDBResidue::residue_atom(Size index) {
        return m_residue_atoms[index];
    };


    // atom finding method
    const PDBAtom& PDBResidue::find_atom(const std::string& atom_name) const {
        auto it = std::find_if(m_residue_atoms.begin(), m_residue_atoms.end(),
                               [&atom_name](const PDBAtom& atom){
                                   return (atom.name() == atom_name);
                               });
        return (*it);
    };


    // residue properties accessors
    Size PDBResidue::n_atoms() const {
        return m_residue_atoms.size();
    };


    // instantiation from a PDBAtomRecord
    PDBResidue PDBResidue::create_from_record(const PDBAtomRecord& pdb_record) {
        PDBResidue pdb_residue;
        pdb_residue.m_residue_serial = pdb_record._residue_seq._record_value;
        pdb_residue.m_residue_name = pdb_record._residue_name._record_value;
        return pdb_residue;
    };


}
