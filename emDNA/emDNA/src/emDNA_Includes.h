// emDNA Includes
// Nicolas Clauvelin


#ifndef emDNA_emDNA_Includes_h
#define emDNA_emDNA_Includes_h


// DNASim headers
// -- maths
#include <Vector3.h>
#include <VectorN.h>
#include <Matrix3.h>
#include <MatrixN.h>
#include <SparseMatrixEntries.h>
#include <SparseMatrixN.h>
#include <Triad.h>
#include <AffineTransformation.h>
#include <Segment.h>
#include <CurveTopology.h>
#include <DiscreteRibbon.h>
#include <BpStepTwistDensity.h>
#include <PointCloudKdTree.h>
// --- DNA
#include <StepParameters.h>
#include <Sequence.h>
#include <StepParametersDB.h>
#include <ForceConstantsDB.h>
// --- I/O
#include <EnhancedString.h>
#include <InputFileHandler.h>
#include <OutputFileHandler.h>
#include <DateParser.h>
// --- utils
#include <OptionsManager.h>
#include <Parser_x3DNA.h>
// --- namespace
using namespace DNASim;


// types renaming
using BasePair = Triad;
using BasePairConstIt = std::vector<BasePair>::const_iterator;
using SizePair = std::pair<Size,Size>;
using EnergyContrib = std::pair<std::string,Real>;


// global constants
namespace emDNAConstants {

    // step parameters dimension
    const Size StepParametersDim = 6;

    // step and force dbs dimension
    const Size SeqDim = 4;

    // pulling force scaling
    const Real PullingForceScaling = Real(1./41.0);

}


// dofs trimming size for the various boundary conditions
namespace BoundaryConditionsDofsTrimming {
    const Size FreeCollection = 0;
    const Size EEDCollection = 3;
    const Size EEDRCollection = 6;
}


// alglib minimizer settings data structure
struct AlglibMinSettings {
    Size _max_iterations;               // maximum number of iterations
    Real _threshold_dx;                 // threshold on relative step size
    Real _threshold_f;                  // threshold on function improvement
    Real _threshold_g;                  // threshold on gradient norm
    Real _max_step_size;                // upper bound for step size
    AlglibMinSettings() :
    _max_iterations(Size(0)),
    _threshold_dx(FLOAT_INIT),
    _threshold_f(FLOAT_INIT),
    _threshold_g(FLOAT_INIT),
    _max_step_size(FLOAT_INIT) {};
};


// callback options
struct CallbackOptions {
    bool _callback_output = false;
    bool _output_current_energy = false;
    Size _output_collection_frequency = 0;
};


// specialized print_vector function for SizePair
namespace DNASim {
    template <>
    inline void print_vector<SizePair>(const std::vector<SizePair>& data,
                                       const char& opening,
                                       const char& closing,
                                       const char& separator) {

        Size n = data.size();
        if (n == 0)
            return;

        std::cout << opening;
        for (Size i=0; i<n; ++i) {
            std::cout << data[i].first << ":" << data[i].second;
            if (i != n-1)
                std::cout << separator;
        };
        std::cout << closing;
        
    };
}


#endif  // emDNA_emDNA_Includes_h
