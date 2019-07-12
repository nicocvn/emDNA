// PDBAtomRecord structure
// Nicolas Clauvelin


// this file declares and implements the data structures to deal with PDB ATOM
// records

// this is ugly ... thumbs up PDB ;)


#ifndef DNASim_PDBAtomRecord_h
#define DNASim_PDBAtomRecord_h


#include <EnhancedString.h>


// PDB format
//
//COLUMNS        DATA  TYPE    FIELD        DEFINITION
//-------------------------------------------------------------------------------------
//1 -  6         Record name   "ATOM  "
//7 - 11         Integer       serial       Atom  serial number.
//13 - 16        Atom          name         Atom name.
//17             Character     altLoc       Alternate location indicator.
//18 - 20        Residue name  resName      Residue name.
//22             Character     chainID      Chain identifier.
//23 - 26        Integer       resSeq       Residue sequence number.
//27             AChar         iCode        Code for insertion of residues.
//31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
//39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
//47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
//55 - 60        Real(6.2)     occupancy    Occupancy.
//61 - 66        Real(6.2)     tempFactor   Temperature  factor.
//77 - 78        LString(2)    element      Element symbol, right-justified.
//79 - 80        LString(2)    charge       Charge  on the atom.


// field widths macros
#define RECORD_NAME_SZ 6
#define ATOM_SERIAL_SZ 5
#define ATOM_NAME_SZ 4
#define ATOM_ALTLOC_SZ 1
#define RES_NAME_SZ 3
#define CHAIN_ID_SZ 1
#define RES_SEQ_SZ 4
#define ICODE_SZ 1
#define ATOM_COORD_SZ 8
#define OCC_SZ 6
#define TEMP_FACTOR_SZ 6

// field start columns
// minus one to account for numbering starting with zero
#define RECORD_NAME_LOC 0
#define ATOM_SERIAL_LOC 6
#define ATOM_NAME_LOC 12
#define ATOM_ALTLOC_LOC 16
#define RES_NAME_LOC 17
#define CHAIN_ID_LOC 21
#define RES_SEQ_LOC 22
#define ICODE_LOC 26
#define ATOM_COORDX_LOC 30
#define ATOM_COORDY_LOC 38
#define ATOM_COORDZ_LOC 46
#define OCC_LOC 54
#define TEMP_FACTOR_LOC 60


namespace DNASim {


    template <class RecordType> struct RecordField {
        Size _field_size;
        RecordType _record_value;
        RecordField(Size field_size) :
        _field_size(field_size),
        _record_value(RecordType()) {};
    };


    struct PDBAtomRecord {

        // atom serial number
        RecordField<Integer> _atom_serial =
        RecordField<Integer>(ATOM_SERIAL_SZ);

        // atom name
        RecordField<std::string> _atom_name =
        RecordField<std::string>(ATOM_NAME_SZ);

        // atom alternate location indicator
        RecordField<std::string> _alternate_loc =
        RecordField<std::string>(ATOM_ALTLOC_SZ);

        // residue name
        RecordField<std::string> _residue_name =
        RecordField<std::string>(RES_NAME_SZ);

        // chain identifier
        RecordField<std::string> _chain_id =
        RecordField<std::string>(CHAIN_ID_SZ);

        // residue sequence number
        RecordField<Integer> _residue_seq = RecordField<Integer>(RES_SEQ_SZ);

        // code for insertion of residues
        RecordField<std::string> _i_code = RecordField<std::string>(ICODE_SZ);

        // atom coordinates Real(8.3)
        RecordField<Real> _atom_X = RecordField<Real>(ATOM_COORD_SZ);
        RecordField<Real> _atom_Y = RecordField<Real>(ATOM_COORD_SZ);
        RecordField<Real> _atom_Z = RecordField<Real>(ATOM_COORD_SZ);

        // occupancy Real(6.2)
        RecordField<Real> _occupancy = RecordField<Real>(OCC_SZ);

        // temperature factor Real(6.2)
        RecordField<Real> _temp_factor = RecordField<Real>(TEMP_FACTOR_SZ);

        // PDB formatted output
        std::string pdb_ATOM_record_string() const {

            // stringstream
            std::stringstream record;

            // record type
            record.width(RECORD_NAME_SZ);
            record << std::left << "ATOM";

            // atom serial
            record.width(ATOM_SERIAL_SZ);
            record << std::right << _atom_serial._record_value;

            // blank space
            record << " ";

            // atom name
            record.width(ATOM_NAME_SZ);
            record << std::left << _atom_name._record_value;

            // atom alternate location
            record.width(ATOM_ALTLOC_SZ);
            record << _alternate_loc._record_value;

            // residue name
            record.width(RES_NAME_SZ);
            record << std::right << _residue_name._record_value;

            // blank space
            record << " ";

            // chain id
            record.width(CHAIN_ID_SZ);
            record << _chain_id._record_value;

            // residue sequence number
            record.width(RES_SEQ_SZ);
            record << std::right << _residue_seq._record_value;

            // icode
            record.width(ICODE_SZ);
            record << _i_code._record_value;

            // blank spaces
            record << "   ";

            // atom coordinates
            record.width(ATOM_COORD_SZ);
            record.setf(std::ios::fixed, std::ios::floatfield);
            record.precision(3);
            record << std::right << _atom_X._record_value;
            record.width(ATOM_COORD_SZ);
            record.setf(std::ios::fixed, std::ios::floatfield);
            record.precision(3);
            record << std::right << _atom_Y._record_value;
            record.width(ATOM_COORD_SZ);
            record.setf(std::ios::fixed, std::ios::floatfield);
            record.precision(3);
            record << std::right << _atom_Z._record_value;

            // occupancy
            record.width(OCC_SZ);
            record.setf(std::ios::fixed, std::ios::floatfield);
            record.precision(2);
            record << std::right << _occupancy._record_value;

            // temperature factor
            record.width(TEMP_FACTOR_SZ);
            record.setf(std::ios::fixed, std::ios::floatfield);
            record.precision(2);
            record << std::right << _temp_factor._record_value;

            return record.str();

        };

        // static method for initialization from a ATOM record string
        static PDBAtomRecord
        init_from_ATOM_record_string(const std::string& record_str) {

            // lambda function for record extraction
            auto extract_record = [&](Size record_loc, Size record_width) {
                return
                EnhancedString::remove_spaces(std::string(record_str,
                                                          record_loc,
                                                          record_width));
            };

            // new PDBAtomRecord
            PDBAtomRecord pdb_atom_record;

            // check for record type
            std::string record_type = extract_record(RECORD_NAME_LOC,
                                                     RECORD_NAME_SZ);
            assert(record_type == "ATOM");

            // atom serial
            std::string atom_serial = extract_record(ATOM_SERIAL_LOC,
                                                     ATOM_SERIAL_SZ);
            pdb_atom_record._atom_serial._record_value =
            EnhancedString::convert_from_string<Integer>(atom_serial);

            // atom name
            std::string atom_name = extract_record(ATOM_NAME_LOC,
                                                   ATOM_NAME_SZ);
            pdb_atom_record._atom_name._record_value = atom_name;

            // alternate location
            std::string atom_altloc = extract_record(ATOM_ALTLOC_LOC,
                                                     ATOM_ALTLOC_SZ);
            pdb_atom_record._alternate_loc._record_value = atom_altloc;

            // residue name
            std::string res_name = extract_record(RES_NAME_LOC,
                                                  RES_NAME_SZ);
            pdb_atom_record._residue_name._record_value = res_name;

            // chain identifier
            std::string chain_id = extract_record(CHAIN_ID_LOC,
                                                  CHAIN_ID_SZ);
            pdb_atom_record._chain_id._record_value = chain_id;

            // residue seq
            std::string res_seq = extract_record(RES_SEQ_LOC,
                                                 RES_SEQ_SZ);
            pdb_atom_record._residue_seq._record_value =
            EnhancedString::convert_from_string<Integer>(res_seq);

            // icode
            std::string icode = extract_record(ICODE_LOC,
                                               ICODE_SZ);
            pdb_atom_record._i_code._record_value = icode;

            // atom coordinates
            std::string atom_x = extract_record(ATOM_COORDX_LOC,
                                                 ATOM_COORD_SZ);
            pdb_atom_record._atom_X._record_value =
            EnhancedString::convert_from_string<Real>(atom_x);
            std::string atom_y = extract_record(ATOM_COORDY_LOC,
                                                ATOM_COORD_SZ);
            pdb_atom_record._atom_Y._record_value =
            EnhancedString::convert_from_string<Real>(atom_y);
            std::string atom_z = extract_record(ATOM_COORDZ_LOC,
                                                ATOM_COORD_SZ);
            pdb_atom_record._atom_Z._record_value =
            EnhancedString::convert_from_string<Real>(atom_z);

            // occupancy
            std::string occupancy = extract_record(OCC_LOC,
                                                   OCC_SZ);
            pdb_atom_record._occupancy._record_value =
            EnhancedString::convert_from_string<Real>(occupancy);

            // temperature factor
            std::string temp_f = extract_record(TEMP_FACTOR_LOC,
                                                TEMP_FACTOR_SZ);
            pdb_atom_record._temp_factor._record_value =
            EnhancedString::convert_from_string<Real>(temp_f);

            return pdb_atom_record;

        };


        PDBAtomRecord() = default;
        ~PDBAtomRecord() = default;
        PDBAtomRecord(const PDBAtomRecord& pdb_atom_record) = default;
        PDBAtomRecord(PDBAtomRecord&& pdb_atom_record) = default;
        PDBAtomRecord&
        operator=(const PDBAtomRecord& pdb_atom_record) = default;
        PDBAtomRecord&
        operator=(PDBAtomRecord&& pdb_atom_record) = default;


    };


}


#undef RECORD_NAME_SZ
#undef ATOM_SERIAL_SZ
#undef ATOM_NAME_SZ
#undef ATOM_ALTLOC_SZ
#undef RES_NAME_SZ
#undef CHAIN_ID_SZ
#undef RES_SEQ_SZ
#undef ICODE_SZ
#undef ATOM_COORD_SZ
#undef OCC_SZ
#undef TEMP_FACTOR_SZ

#undef RECORD_NAME_LOC
#undef ATOM_SERIAL_LOC
#undef ATOM_NAME_LOC
#undef ATOM_ALTLOC_LOC
#undef RES_NAME_LOC
#undef CHAIN_ID_LOC
#undef RES_SEQ_LOC
#undef ICODE_LOC
#undef ATOM_COORDX_LOC
#undef ATOM_COORDY_LOC
#undef ATOM_COORDZ_LOC
#undef OCC_LOC
#undef TEMP_FACTOR_LOC


#endif  // DNASim_PDBAtomRecord_h
