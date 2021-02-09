// BpCollection_Interface class
// Nicolas Clauvelin


// abstract interface for the computation of energy and gradient of a bp
// collection

// this class needs to be derived depending on the type of bp collection:
// with or w/o frozen steps, anchored, clamped, ...

// electrostatics interactions are taken into account or not by using the
// toggle_on/off_electrostatics() methods


#ifndef emDNA_BpCollection_Interface_h
#define emDNA_BpCollection_Interface_h


#include "DNA/IndexManager.h"
#include "DNA/BpCollection.h"


class BpCollection_Interface {


public:

    // bp collection accessor/modifier
    const BpCollection& bp_collection() const;
    virtual void set_bp_collection(const BpCollection& bp_collection);

    // electrostatic interactions flag modifiers
    void toggle_on_electrostatics();
    void toggle_off_electrostatics();

    // bp collection energy
    // this is the objective function
    virtual Real bp_collection_energy() const;

    // bp collection energy contribs
    // returns the list of contributions to the energy
    virtual std::vector<EnergyContrib> energy_contribs() const;

    // bp collection free dofs accessor
    // this is the set of dofs on which the minimization is performed
    // this method returns all free dofs independently of the boundary
    // conditions (needs to be override in derived classes)
    virtual VectorN bp_collection_inline_free_dofs() const;

    // bp collection free dofs scalings
    // this is the scalings for the free dofs used in the minimization
    virtual VectorN bp_collection_free_dofs_scalings() const;

    // bp collection elastic energy gradient
    // this is the gradient of the elastic energy only
    VectorN bp_collection_elastic_energy_gradient() const;

    // bp collection free dofs gradient
    // this is the gradient used for the minimization
    virtual VectorN bp_collection_free_dofs_gradient() const =0;

    // update method
    virtual void update_with_new_inline_free_dofs(const VectorN&
                                                  inline_free_dofs) =0;

    // virtual destructor
    virtual ~BpCollection_Interface() {};


protected:

    // constructors
    BpCollection_Interface();
    BpCollection_Interface(const BpCollection_Interface&
                           bp_collection_interface);

    // copy operator
    BpCollection_Interface& operator=(const BpCollection_Interface&
                                      bp_collection_interface);

    // method for rebuilding dofs from free dofs
    // dofs for the step in the range [0, end) are rebuild
    std::vector<BpStepDofs>
    rebuild_dofs_from_free_dofs(const Size& range_end,
                                const std::vector<BpStepDofs>& free_dofs) const;

    IndexManager m_idx_mgr;             // index manager
    BpCollection m_bp_collection;       // bp collection

    // boolean flag for taking into account electrostatic interactions
    // true -> electrostatic included
    // false -> electrostatic not included
    // by default we set it to false
    bool m_elec_flag;


};


// static repacking method
// for incomplete inline values the corresponding dofs are set to zero
// incomplete means not a multiple of emDNAConstants::StepParametersDim
template<class T>
static std::vector<T> repack_inline_vector(const VectorN& values,
                                           Size pack_size =
                                           emDNAConstants::StepParametersDim) {

    // container
    std::vector<T> packed;

    // size info
    div_t modulo = div((Integer)values.size(), pack_size);
    Size n_full_steps = modulo.quot;
    packed.reserve(n_full_steps);

    // loop for full steps
    for (Size i=0; i<n_full_steps; ++i)
        packed.push_back(T( VectorN(values.slice(pack_size*i,
                                                 pack_size*(i+1))) ));

    // remaining values
    if (modulo.rem != 0) {
        VectorN final_slice(pack_size, FLOAT_INIT);
        final_slice.set_slice(0,
                              values.slice(pack_size*n_full_steps,
                                           values.size()));
        packed.push_back(T(final_slice));
    };

    return packed;

};


#endif  // emDNA_BpCollection_Interface_h
