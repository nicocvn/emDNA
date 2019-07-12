// FreeBpCollection class
// Nicolas Clauvelin


// bp collection interface for free bp collection:
//  - no end-to-end conditions

// the minimization of such a collection with respect to the dofs is
// somehow ill-defined (more efficient to minimize with respect to the step
// parameters directly)
// therefore, this implementation is for testing purposes only


#ifndef emDNA_FreeBpCollection_h
#define emDNA_FreeBpCollection_h


#include <BpCollection_Interface.h>


class FreeBpCollection final : public BpCollection_Interface {


public:

    // constructors
    FreeBpCollection();
    FreeBpCollection(const FreeBpCollection& free_collection);
    ~FreeBpCollection();

    // copy operator
    FreeBpCollection& operator=(const FreeBpCollection&
                                free_collection);

    // bp collection free dofs gradient
    // this is the gradient used for the minimization
    VectorN bp_collection_free_dofs_gradient() const override;

    // update method
    void
    update_with_new_inline_free_dofs(const VectorN& inline_free_dofs) override;


private:

    // frozen steps rebuilding method
    std::vector<BpStepDofs>
    rebuild_full_dofs_from_free_dofs(const std::vector<BpStepDofs>&
                                     free_dofs) const;


};


#endif  // emDNA_FreeBpCollection_h
