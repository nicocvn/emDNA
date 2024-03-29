// ClampedBpCollection class
// Nicolas Clauvelin


// bp collection interface for clamped bp collection:
//  - imposed end-to-end distance and rotation


#ifndef emDNA_ClampedBpCollection_h
#define emDNA_ClampedBpCollection_h


#include "DNA/BpCollection_Interface.h"


class ClampedBpCollection final : public BpCollection_Interface {


public:

    // constructors
    ClampedBpCollection();
    ClampedBpCollection(const ClampedBpCollection& clamped_collection);
    ~ClampedBpCollection();

    // copy operator
    ClampedBpCollection& operator=(const ClampedBpCollection&
                                   clamped_collection);

    // bp collection free dofs accessor
    // this is the set of dofs on which the minimization is performed
    // for EEDR boundary conditions this is all free dofs minus the dofs of the
    // bc step
    VectorN bp_collection_inline_free_dofs() const override;

    // bp collection free dofs scalings
    // this is the scalings for the free dofs used in the minimization
    VectorN bp_collection_free_dofs_scalings() const override;

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


#endif  // emDNA_ClampedBpCollection_h
