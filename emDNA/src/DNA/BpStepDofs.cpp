// BpStepDofs class
// Nicolas Clauvelin


#include <BpStepDofs.h>
#include <BpGeometryFunctions.h>


// class default constructor
BpStepDofs::BpStepDofs() :
m_dofs(emDNAConstants::StepParametersDim, FLOAT_INIT) {};


// class constructor with initialization from two frames
BpStepDofs::BpStepDofs(const BasePair& bp1, const BasePair& bp2) {
    *this = BpStepDofs::bp_step_dofs(bp1, bp2);
};


// class constructor with initialization from a vector
BpStepDofs::BpStepDofs(const VectorN& dofs_values) : m_dofs(dofs_values) {

    // size checking
    DS_ASSERT(m_dofs.size()==emDNAConstants::StepParametersDim,
              "initializing dofs with a ill-sized vector\n"
              "vector size is " << dofs_values.size());

};


// class constructor with initialization from a string
BpStepDofs::BpStepDofs(const std::string& string_dofs) : m_dofs(string_dofs) {

    // size checking
    DS_ASSERT(m_dofs.size()==emDNAConstants::StepParametersDim,
              "initializing dofs with a ill-sized vector string\n"
              "vector string is " << string_dofs);

};


// base pair rebuild method
BasePair BpStepDofs::rebuild_step_last_base_pair(const BasePair& bp) const {

    // ZYZ Euler angles
    const BpGeometryFunctions::ThetaAngles theta_angles(m_dofs[0]*RAD_2_DEG,
                                                        m_dofs[1]*RAD_2_DEG,
                                                        m_dofs[2]*RAD_2_DEG);
    const BpGeometryFunctions::ZYZEulerAngles
    euler_angles(BpGeometryFunctions::Theta_2_ZYZEuler(theta_angles));
    const Vector3 euler_angles_vec(euler_angles._ZetaRadian,
                                   euler_angles._KappaRadian,
                                   euler_angles._EtaRadian);

    // translation
    const Vector3 r(m_dofs[3], m_dofs[4], m_dofs[5]);

    // step last base pair
    BasePair last_bp(bp);
    last_bp.transform(Triad::euler_zyz_transformation(bp, euler_angles_vec));
    last_bp.set_origin(bp.origin()+r);

    // orthogonalize
    last_bp.orthogonalize();

    return last_bp;

};


// static  computation method
BpStepDofs BpStepDofs::bp_step_dofs(const BasePair& bp1, const BasePair& bp2) {

    // theta angles
    const Vector3 euler_angles_vec(Triad::euler_zyz_angles(bp1, bp2));
    BpGeometryFunctions::ZYZEulerAngles euler_angles(euler_angles_vec[X],
                                                     euler_angles_vec[Y],
                                                     euler_angles_vec[Z]);
    const BpGeometryFunctions::ThetaAngles
    theta_angles(BpGeometryFunctions::ZYZEuler_2_Theta(euler_angles));

    // step dofs
    BpStepDofs dofs;
    dofs.value(TILTrad) = theta_angles._Theta1Degree*DEG_2_RAD;
    dofs.value(ROLLrad) = theta_angles._Theta2Degree*DEG_2_RAD;
    dofs.value(TWISTrad) = theta_angles._Theta3Degree*DEG_2_RAD;
    dofs.value(R1) = (bp2.origin()-bp1.origin())[X];
    dofs.value(R2) = (bp2.origin()-bp1.origin())[Y];
    dofs.value(R3) = (bp2.origin()-bp1.origin())[Z];

    return dofs;

};
std::vector<BasePair>
BpStepDofs::rebuild_bps(const std::vector<BpStepDofs>& dofs,
                        const BasePair& first_bp) {

    // number of base pairs
    const Size n_bp = dofs.size()+1;

    // allocation and first base pair
    std::vector<BasePair> base_pairs;
    base_pairs.reserve(n_bp);
    base_pairs.push_back(first_bp);

    // rebuild
    BpStepDofsConstIt step_end = dofs.end();
    for (auto  step_it = dofs.begin(); step_it != step_end; ++step_it) {
        base_pairs.push_back(step_it->
                             rebuild_step_last_base_pair(base_pairs.back()));
    };

    return base_pairs;

};
std::vector<BpStepDofs>
BpStepDofs::compute_dofs(const std::vector<BasePair>& bps) {

    // number of steps
    const Size n_step = bps.size()-1;

    // allocation
    std::vector<BpStepDofs> dofs;
    dofs.reserve(n_step);

    // computation
    BasePairConstIt bp_end = bps.end();
    for (auto bp_it = bps.begin(); bp_it != bp_end-1; ++bp_it)
        dofs.push_back(BpStepDofs::bp_step_dofs(*(bp_it), *(bp_it+1)));

    return dofs;

};
