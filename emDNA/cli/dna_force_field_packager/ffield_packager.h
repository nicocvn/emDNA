// ffield_packager tool
// Nicolas Clauvelin


// this tool is used to package a force field into a binary format (ff
// extension)
// a force field is described by two files:
// - a file listing the intrinsic step parameters,
// - a file describing the force constants matrices.
//
// those files have to contain the data for all the step sequences; refer
// to StepParameters_IdealDNA.h and ForceConstants_IdealDNA.h for more details


#ifndef DNASim_ffield_packager_h
#define DNASim_ffield_packager_h


#include <DNASim.h>
using namespace DNASim;


// data structure
struct CLData {
    std::string _model_name;
    std::string _steps_filename;
    std::string _fmat_filename;
};


// types renaming
template <class T>
using TupleT = std::tuple<BaseSymbol, BaseSymbol, BaseSymbol, BaseSymbol, T>; //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
template <class T>
using TupleTVector = std::vector<TupleT<T>>;


// entry point
int main(int argc, char* argv[]);


// command line parsing option
CLData parse_command_line(int argc, char* argv[]);


// building functions
TupleTVector<VectorN> parse_text_file(const std::string& filename);
StepParametersDB build_steps_db(const std::string& model_name,
                                const TupleTVector<VectorN>& data);
ForceConstantsDB build_fmat_db(const std::string& model_name,
                               const TupleTVector<VectorN>& data);


// force field archive function
template <class ArchiveFormat>
void archive_force_field(const std::string& filename,
                         const StepParametersDB& step_db,
                         const ForceConstantsDB& fmat_db) {

    // output archive
    OutputArchive<ArchiveFormat> archive(filename, false);

    // force field packing
    archive.save(step_db);
    archive.save(fmat_db);
    
};


#endif  // DNASim_ffield_packager_h
