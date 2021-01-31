// PDBParser class
// Nicolas Clauvelin


#include "file_io/InputFileHandler.h"
#include "pdb/PDBAtomRecord.h"
#include "pdb/PDBParser.h"


namespace DNASim {


    // static method for the instantiation of a PDBStructure from a PDB file
    PDBStructure PDBParser::
    create_PDB_structure_from_file(const std::string& pdb_file_name) {

        // file slurping
        std::vector<std::string> pdb_content;
        pdb_content = InputFileHandler::read_file(pdb_file_name, ".");

        // keep only ATOM records
        std::vector<std::string> pdb_atom_record_lines;
        std::for_each(pdb_content.begin(), pdb_content.end(),
                      [&pdb_atom_record_lines](const std::string& line){
#define ATOM_RECORD_LENGTH 4
                          std::string record_type(line, 0, ATOM_RECORD_LENGTH);
#undef ATOM_RECORD_LENGTH
                          if (record_type == "ATOM")
                              pdb_atom_record_lines.push_back(line);
                      });

        // create atom record structures
        PDBAtomRecordVector pdb_atom_records;
        pdb_atom_records.reserve(pdb_atom_record_lines.size());
        std::for_each(pdb_atom_record_lines.begin(),
                      pdb_atom_record_lines.end(),
                      [&pdb_atom_records](const std::string& line){
                          pdb_atom_records.
                          push_back(PDBAtomRecord::
                                    init_from_ATOM_record_string(line));
                      });

        // list of chain identifiers
        std::vector<std::string> all_chain_ids = chain_ids(pdb_atom_records);

        // chain creation
        std::vector<PDBChain> pdb_chains;
        for (const std::string& chain_id : all_chain_ids)
            pdb_chains.push_back(make_chain(chain_id, pdb_atom_records));

        // structure final assembly
        PDBStructure pdb_structure;
        std::for_each(pdb_chains.begin(), pdb_chains.end(),
                      [&pdb_structure](const PDBChain& chain){
                          pdb_structure.add_chain(chain);
                      });

        return pdb_structure;

    };


    // get list of chain identifiers from a list of records
    std::vector<std::string>
    PDBParser::chain_ids(const PDBAtomRecordVector& pdb_records) {

        // list of all chains
        std::vector<std::string> all_chain_ids;
        all_chain_ids.reserve(pdb_records.size());
        std::for_each(pdb_records.begin(), pdb_records.end(),
                      [&all_chain_ids](const PDBAtomRecord& record) {
                          all_chain_ids.push_back(record.
                                                  _chain_id._record_value);
                      });

        // sort and remove duplicates
        // we should not need that ... but PDB format is weird sometimes
        std::sort(all_chain_ids.begin(),all_chain_ids.end());
        all_chain_ids.erase(std::unique(all_chain_ids.begin(),
                                        all_chain_ids.end()),
                            all_chain_ids.end());

        return all_chain_ids;
        
    };


    // create a chain from a list of records for a specific chain identifier
    PDBChain PDBParser::make_chain(const std::string& chain_id,
                                   const PDBAtomRecordVector& pdb_records) {

        // residue container
        PDBChain pdb_chain;
        pdb_chain.set_identifier(chain_id);

        // lambda function for testing the presence of a residue in the
        // container
        auto find_residue_in_chain =
        [&pdb_chain](const PDBAtomRecord& record) -> bool {
            auto it =
            std::find_if(pdb_chain.residue_begin(), pdb_chain.residue_end(),
                         [&record](const PDBResidue& res){
                             return (res.residue_sequence_index() ==
                                     record._residue_seq._record_value);
                         });
            return (it == pdb_chain.residue_end()) ? false : true;
        };

        // loop over the records to create the residues
        std::for_each(pdb_records.begin(), pdb_records.end(),
                      [&chain_id, &pdb_chain, &find_residue_in_chain]
                      (const PDBAtomRecord& record) {

                          // check for the chain
                          if (chain_id == record._chain_id._record_value) {

                              // check if the residue has already been created
                              // if not we create it
                              if (!find_residue_in_chain(record)) {
                                  pdb_chain.
                                  add_residue(PDBResidue::
                                              create_from_record(record));
                              };

                              // sequence index
                              Size seq_idx = record._residue_seq._record_value;

                              // add the atom to the chain
                              pdb_chain.
                              chain_residue_from_sequence_index(seq_idx).
                              add_atom(PDBAtom::create_from_record(record));

                          };

                      });

        return pdb_chain;

    };


}
