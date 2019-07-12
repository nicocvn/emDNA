// IndexManager class
// Nicolas Clauvelin


// index manager implementation for the bp collection

// this class is used to manipulate the indexes of the bp step
// free steps = every non-frozen steps
// frozen steps = steps for which the step parameters are imposed


#ifndef emDNA_IndexManager_h
#define emDNA_IndexManager_h


#include <BpStepDofs.h>
#include <BpStepParams.h>
#include <IndexPair.h>


// iterators typedefs
typedef std::vector<IndexPair>::const_iterator IndexPairConstIt;


class IndexManager {


public:

    // constructor
    IndexManager();
    IndexManager(const IndexManager& idx_mgr);
    ~IndexManager();

    // copy operator
    IndexManager& operator=(const IndexManager& idx_mgr);

    // setup methods
    // the frozen domains are given with inclusive intervals [i,j]
    void set_n_of_bp_step(Size n_bp_steps);
    void set_frozen_steps_indexes(const std::vector<SizePair>&
                                  frozen_step_domains);

    // size accessors
    Size n_of_dofs_variables() const;
    Size n_of_free_dofs_variables() const;
    Size n_of_parameters_variables() const;
    Size n_of_free_parameters_variables() const;

    // free step indexes accessors
    Size n_of_free_steps() const;
    std::vector<Size> inline_free_dofs_indices() const;
    std::vector<Size> free_step_collection_indices() const;
    IndexPairConstIt free_step_indexes_begin() const;
    IndexPairConstIt free_step_indexes_end() const;


    // bc reduced step index accessor
    // this corresponds to the last free step
    const IndexPair& bc_reduced_step() const;

    // frozen step indexes accessors
    Size n_of_frozen_steps() const;
    std::vector<Size> frozen_step_collection_indices() const;
    IndexPairConstIt frozen_step_indexes_begin() const;
    IndexPairConstIt frozen_step_indexes_end() const;

    // dof and parameter global coordinate methods
    Size dof_global_coordinate(Size step_index,
                               BpStepDof dof_index) const;
    Size parameter_global_coordinate(Size step_index,
                                     BpStepParameter param_index) const;


private:

    // total number of steps
    Size m_n_steps;

    // free and frozen step indexes pairs
    std::vector<IndexPair> m_free_steps;
    std::vector<IndexPair> m_frozen_steps;


};


#endif  // emDNA_IndexManager_h
