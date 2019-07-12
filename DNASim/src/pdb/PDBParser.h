// PDBParser class
// Nicolas Clauvelin


// this class implements a basic PDB parser (no guarantee)

// limitations:
// - no model implementations
// - parse ATOM records only


#ifndef DNASim_PDBParser_h
#define DNASim_PDBParser_h


#include <PDBStructure.h>


namespace DNASim {


    class PDBAtomRecord;
    typedef std::vector<PDBAtomRecord> PDBAtomRecordVector;


    class PDBParser {


    public:

        // static method for the instantiation of a PDBStructure from a PDB file
        static PDBStructure
        create_PDB_structure_from_file(const std::string& pdb_file_name);


    private:

        // private parsing methods
        static std::vector<std::string> chain_ids(const PDBAtomRecordVector&
                                                  pdb_records);
        static PDBChain make_chain(const std::string& chain_id,
                                   const PDBAtomRecordVector& pdb_records);

        // deleted constructors and operators
        PDBParser() = delete;
        ~PDBParser() = delete;
        PDBParser(const PDBParser& pdb_parser) = delete;
        PDBParser(PDBParser&& pdb_parser) = delete;
        PDBParser& operator=(const PDBParser& pdb_parser) = delete;
        PDBParser operator=(PDBParser&& pdb_parser) = delete;


    };


}


#endif  // DNASim_PDBParser_h
