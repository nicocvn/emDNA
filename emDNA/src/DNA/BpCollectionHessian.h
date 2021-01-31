// BpCollectionHessian functions
// Nicolas Clauvelin


// implementation of Hessian computations

#ifndef emDNA_BpCollectionHessian_h
#define emDNA_BpCollectionHessian_h


#include "DNA/StepBlock.h"
class BpCollection;
class IndexManager;


namespace BpCollectionHessian {


    // free collection Hessian matrix
    MatrixN free_collection_Hessian(const BpCollection& bp_collection,
                                    const IndexManager& idx_mgr);

//    // following code is in an unnamed namespace
//    namespace {


    // EEDR boundary conditions B matrix as StepBlock
    StepBlock eedr_B_matrix(const Size& step_index,
                            const BpCollection& bp_collection,
                            const IndexManager& idx_mgr);


    // first order eedr contributions
    void order1_eerd_contributions(StepBlockArray2D& free_hessian,
                                   const BpCollection& bp_collection,
                                   const IndexManager& idx_mgr);


//    }


}


#endif  // emDNA_BpCollectionHessian_h
