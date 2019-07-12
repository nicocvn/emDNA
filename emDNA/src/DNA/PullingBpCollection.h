// PullingBpCollection class
// Nicolas Clauvelin


// bp collection interface for pulling bp collection:
//  - no end-to-end rotation
//  - external force applied to the terminal bp


#ifndef emDNA_PullingBpCollection_h
#define emDNA_PullingBpCollection_h


#include <BpCollection_Interface.h>


class PullingBpCollection final : public BpCollection_Interface {


public:

    // constructors
    PullingBpCollection();
    PullingBpCollection(const PullingBpCollection& pulled_collection);
    ~PullingBpCollection();

    // copy operator
    PullingBpCollection& operator=(const PullingBpCollection&
                                   pulled_collection);

    // bp collection modifier
    void set_bp_collection(const BpCollection& bp_collection) override;

    // bp collection energy
    // this is the objective function
    // the energy is always calculated on the set of free steps
    Real bp_collection_energy() const override;

    // bp collection energy contribs
    // returns the list of contributions to the energy
    std::vector<EnergyContrib> energy_contribs() const override;

    // bp collection free dofs gradient
    // this is the gradient used for the minimization
    VectorN bp_collection_free_dofs_gradient() const override;

    // update method
    void
    update_with_new_inline_free_dofs(const VectorN& inline_free_dofs) override;

    // pulling force accessor/modifier
    const Vector3& pulling_force() const;
    void set_pulling_force(const Vector3& f);


private:

    // pulling energy method
    Real pulling_energy() const;

    // pulling force
    Vector3 m_pulling_force;
    Vector3 m_rest_state_endpoint;

    // frozen steps rebuilding method
    std::vector<BpStepDofs>
    rebuild_full_dofs_from_free_dofs(const std::vector<BpStepDofs>&
                                     free_dofs) const;
    
    
};


#endif  // emDNA_PullingBpCollection_h
