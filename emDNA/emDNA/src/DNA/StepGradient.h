// StepGradient structure
// Nicolas Clauvelin


#ifndef emDNA_StepGradient_h
#define emDNA_StepGradient_h


#include <emDNA_Includes.h>
class IndexManager;


// step gradient structure
// the gradient is stored as an angular and translational parts
struct StepGradient {
    Vector3 _rotation;
    Vector3 _translation;
    StepGradient() = default;
    ~StepGradient() = default;
    StepGradient(const StepGradient& sg) = default;
    StepGradient(StepGradient&& sg) = default;
    StepGradient& operator=(const StepGradient& sg) = default;
    StepGradient& operator=(StepGradient&& sg) = default;
    friend std::ostream& operator<<(std::ostream& os,
                                    const StepGradient& g) {
        os << "rotation = " << g._rotation;
        os << "\t";
        os << "translation = " << g._translation;
        return os;
    };
};
using StepGradientVec = std::vector<StepGradient>;


// flatenning function for a vector of StepGradient
VectorN flatten_StepGradientVec(const StepGradientVec& v);


// filtering function for a vector of StepGradient
StepGradientVec filter_free_dofs_grads(const StepGradientVec& v,
                                       const IndexManager& idx_mgr);


#endif  // emDNA_StepGradient_h
