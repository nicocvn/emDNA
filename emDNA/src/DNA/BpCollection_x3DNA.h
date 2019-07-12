// BpCollection_x3DNA class
// Nicolas Clauvelin


#ifndef emDNA_BpCollection_x3DNA_h
#define emDNA_BpCollection_x3DNA_h


#include <emDNA_Includes.h>
class BpCollection;


class BpCollection_x3DNA {


public:

    // x3DNA reading methods
    static BpCollection
    create_bp_collection_from_step_parameters_file(const std::string& filename);
    static BpCollection
    create_bp_collection_from_base_pairs_file(const std::string& filename);

    // x3DNA output methods
    static void
    format_as_step_parameters_file(std::stringstream& output,
                                   const BpCollection& bp_collection);
    static void
    format_as_base_pairs_file(std::stringstream& output,
                              const BpCollection& bp_collection);

    // x3DNA writing methods
    static void
    write_as_step_parameters_file(const std::string& filename,
                                  const BpCollection& bp_collection);
    static void
    write_as_base_pairs_file(const std::string& filename,
                             const BpCollection& bp_collection);


private:

    // constructors and operators
    BpCollection_x3DNA() = delete;
    BpCollection_x3DNA(const BpCollection_x3DNA& bp_x3DNA) = delete;
    BpCollection_x3DNA(BpCollection_x3DNA&& bp_x3DNA) = delete;
    BpCollection_x3DNA& operator=(const BpCollection_x3DNA& bp_x3DNA) = delete;
    BpCollection_x3DNA& operator=(BpCollection_x3DNA&& bp_x3DNA) = delete;
    ~BpCollection_x3DNA() = delete;


};


#endif  // emDNA_BpCollection_x3DNA_h
