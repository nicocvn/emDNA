// BpStepDofs class
// Nicolas Clauvelin


// base-pair step degrees of freedom (dofs) class

// base-pair step dofs units:
//  TILT    radian
//  ROLL    radian
//  TWIST   radian
//  R1      angstroms
//  R2      angstroms
//  R3      angstroms


#ifndef emDNA_BpStepDofs_h
#define emDNA_BpStepDofs_h


#include <emDNA_Includes.h>


// enum type for step parameters
enum BpStepDof {
    TILTrad=0, ROLLrad=1, TWISTrad=2,
    R1=3, R2=4, R3=5
};


class BpStepDofs {


public:

    // constructors
    BpStepDofs();
    BpStepDofs(const BpStepDofs& bp_step_dofs) = default;
    BpStepDofs(BpStepDofs&& bp_step_dofs) = default;
    BpStepDofs(const BasePair& bp1, const BasePair& bp2);
    BpStepDofs(const VectorN& dofs_values);
    BpStepDofs(const std::string& string_dofs);
    ~BpStepDofs() = default;

    // copy and move operators
    BpStepDofs& operator=(const BpStepDofs& bp_step_dofs) = default;
    BpStepDofs& operator=(BpStepDofs&& bp_step_dofs) = default;

    // parameters accessors/modifiers
    inline const Real& value(const BpStepDof& i) const {
        return m_dofs[i];
    };
    inline Real& value(const BpStepDof& i) {
        return m_dofs[i];
    };

    // inline vector
    inline const VectorN& inline_vector() const {
        return m_dofs;
    };
    inline VectorN& inline_vector() {
        return m_dofs;
    };

    // base pair rebuild method
    BasePair rebuild_step_last_base_pair(const BasePair& bp) const;

    // static computation methods
    static
    BpStepDofs bp_step_dofs(const BasePair& bp1, const BasePair& bp2);
    static
    std::vector<BasePair> rebuild_bps(const std::vector<BpStepDofs>& dofs,
                                      const BasePair& first_bp);
    static
    std::vector<BpStepDofs> compute_dofs(const std::vector<BasePair>& bps);


private:

    VectorN m_dofs;

    
};


// const iterator typedef
typedef std::vector<BpStepDofs>::const_iterator BpStepDofsConstIt;


#endif  // emDNA_BpStepDofs_h
