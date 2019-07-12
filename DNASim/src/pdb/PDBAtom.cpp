// PDBAtom class
// Nicolas Clauvelin


#include <PDBAtomRecord.h>
#include <PDBAtom.h>


namespace DNASim {


    // atom serial accessor/modifier
    const Integer& PDBAtom::atom_serial() const {
        return m_atom_serial;
    };
    void PDBAtom::set_atom_serial(const Integer& serial) {
        m_atom_serial = serial;
    };


    // atom name accessor/modifier
    const std::string& PDBAtom::name() const {
        return m_atom_name;
    };
    void PDBAtom::set_name(const std::string& atom_name) {
        m_atom_name = atom_name;
    };


    // atom position accessor/modifier
    const Vector3& PDBAtom::position() const {
        return m_atom_position;
    };
    void PDBAtom::set_position(const Vector3& position) {
        m_atom_position = position;
    };


    // instantiation from a PDBAtomRecord
    PDBAtom PDBAtom::create_from_record(const PDBAtomRecord& pdb_record) {
        PDBAtom pdb_atom;
        pdb_atom.m_atom_serial = pdb_record._atom_serial._record_value;
        pdb_atom.m_atom_name = pdb_record._atom_name._record_value;
        pdb_atom.m_atom_position = Vector3(pdb_record._atom_X._record_value,
                                           pdb_record._atom_Y._record_value,
                                           pdb_record._atom_Z._record_value);
        return pdb_atom;
    };


}
