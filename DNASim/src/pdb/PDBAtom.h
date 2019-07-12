// PDBAtom class
// Nicolas Clauvelin


#ifndef DNASim_PDBAtom_h
#define DNASim_PDBAtom_h


#include <Vector3.h>


namespace DNASim {


    class PDBAtomRecord;


    class PDBAtom {
        

    public:

        // constructors
        PDBAtom() = default;
        PDBAtom(const PDBAtom& pdb_atom) = default;
        PDBAtom(PDBAtom&& pdb_atom) = default;
        ~PDBAtom() = default;

        // copy and move operators
        PDBAtom& operator=(const PDBAtom& pdb_atom) = default;
        PDBAtom& operator=(PDBAtom&& pdb_atom) = default;

        // atom serial accessor/modifier
        const Integer& atom_serial() const;
        void set_atom_serial(const Integer& serial);

        // atom name accessor/modifier
        const std::string& name() const;
        void set_name(const std::string& atom_name);

        // atom position accessor/modifier
        const Vector3& position() const;
        void set_position(const Vector3& position);

        // instantiation from a PDBAtomRecord
        static PDBAtom create_from_record(const PDBAtomRecord& pdb_record);


    private:

        Integer m_atom_serial;
        std::string m_atom_name;
        Vector3 m_atom_position;


    };


}


#endif // DNASim_PDBAtom_h
