// emDNA_Serialization class
// Nicolas Clauvelin


// cereal serialization load/save functions for emDNA


#ifndef emDNA_emDNA_Serialization_h
#define emDNA_emDNA_Serialization_h


#include <emDNA_Includes.h>
// --- serialization
#include <Cereal_Includes.h>
#include <DNASimSerialization.h>
// --- data types
#include <BpStepDofs.h>
#include <BpStepParams.h>


// BpStepDofs
template <class ArchiveType>
void save(ArchiveType& archive, const BpStepDofs& dofs) {
    archive(cereal::make_nvp("inline_values", dofs.inline_vector()));
};
template <class ArchiveType>
void load(ArchiveType& archive, BpStepDofs& dofs) {
    VectorN inline_values;
    archive(inline_values);
    dofs = BpStepDofs(inline_values);
};


// BpStepParams
template <class ArchiveType>
void save(ArchiveType& archive, const BpStepParams& p) {
    archive(cereal::make_nvp("inline_values", p.inline_vector()));
};
template <class ArchiveType>
void load(ArchiveType& archive, BpStepParams& p) {
    VectorN inline_values;
    archive(inline_values);
    p = BpStepParams(inline_values);
};


#endif  // emDNA_emDNA_Serialization_h
