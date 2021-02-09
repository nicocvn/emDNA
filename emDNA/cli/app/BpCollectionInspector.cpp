//// BpCollectionInspector class
//// Nicolas Clauvelin
//
//
//#include <BpCollectionGeometry.h>
//#include <BpCollectionGradient.h>
//#include <BpCollectionElasticEnergy.h>
//#include <BpCollectionInspector.h>
//
//
//// static output to file method
//#define N_DIGITS 6
//#define SPACING 15
//void BpCollectionInspector::write_inspector_data(const std::string&
//                                                 output_file,
//                                                 const BpCollection&
//                                                 bp_collection) {
//
//    // data
//    BpCollectionData data = inspect_bp_collection(bp_collection);
//    const Size n_bp = bp_collection.n_of_base_pairs();
//
//    // data organization
//    std::vector< std::vector<std::string> > data_str;
//    for (Size i=0; i<n_bp-1; ++i) {
//        std::vector<std::string> v;
//        v.push_back(EnhancedString::convert_to_string(i));
//        v.push_back(bp_collection.bp_step_sequence(i).str());
//        v.push_back(Parser_x3DNA::
//                    fixed_width_convert_to_string(data._energy_per_step[i],
//                                                  N_DIGITS));
//        v.push_back(Parser_x3DNA::
//                    fixed_width_convert_to_string(data._bend_angle_per_step[i],
//                                                  N_DIGITS));
//        v.push_back(Parser_x3DNA::
//                    fixed_width_convert_to_string(data._twist_sc_per_step[i],
//                                                  N_DIGITS));
////        v.push_back(fixed_width_convert_to_string(data._force_per_bp[i][X],
////                                                  N_DIGITS));
////        v.push_back(fixed_width_convert_to_string(data._force_per_bp[i][Y],
////                                                  N_DIGITS));
////        v.push_back(fixed_width_convert_to_string(data._force_per_bp[i][Z],
////                                                  N_DIGITS));
////        v.push_back(fixed_width_convert_to_string(data._moment_per_bp[i][X],
////                                                  N_DIGITS));
////        v.push_back(fixed_width_convert_to_string(data._moment_per_bp[i][Y],
////                                                  N_DIGITS));
////        v.push_back(fixed_width_convert_to_string(data._moment_per_bp[i][Z],
////                                                  N_DIGITS));
//        data_str.push_back(v);
//    };
//
////    // add data for the last base pair
////    std::vector<std::string> v;
////    v.push_back(EnhancedString::convert_to_string(n_bp-1));
////    v.push_back("--");
////    v.push_back("--");
////    v.push_back("--");
////    v.push_back(fixed_width_convert_to_string(data._force_per_bp[n_bp-1][X],
////                                              N_DIGITS));
////    v.push_back(fixed_width_convert_to_string(data._force_per_bp[n_bp-1][Y],
////                                              N_DIGITS));
////    v.push_back(fixed_width_convert_to_string(data._force_per_bp[n_bp-1][Z],
////                                              N_DIGITS));
////    v.push_back(fixed_width_convert_to_string(data._moment_per_bp[n_bp-1][X],
////                                              N_DIGITS));
////    v.push_back(fixed_width_convert_to_string(data._moment_per_bp[n_bp-1][Y],
////                                              N_DIGITS));
////    v.push_back(fixed_width_convert_to_string(data._moment_per_bp[n_bp-1][Z],
////                                              N_DIGITS));
////    data_str.push_back(v);
//
//    // header
//    std::vector<std::string> header;
//    header.push_back("#INDEX");
//    header.push_back("SEQ");
//    header.push_back("ENERGY(kBT)");
//    header.push_back("BEND(deg)");
//    header.push_back("TWIST_SC(deg)");
////    header.push_back("F1(pN)");
////    header.push_back("F2(pN)");
////    header.push_back("F3(pN)");
////    header.push_back("M1(pN.A)");
////    header.push_back("M2(pN.A)");
////    header.push_back("M3(pN.A)");
//
//    // output stream
//    std::stringstream output;
//
//    // header
//    Parser_x3DNA::fixed_width_vector_formatting(output, header, N_DIGITS, SPACING);
//
//    // stats values
//    for (Size i=0; i<data_str.size(); ++i)
//        Parser_x3DNA::
//        fixed_width_vector_formatting(output, data_str[i], N_DIGITS, SPACING);
//
//    // stats file
//    OutputFileHandler output_f(output_file, ".");
//    output_f.open();
//    output_f.write_stream(output);
//    output_f.close();
//
//};
//#undef N_DIGITS
//#undef SPACING
//
//
//// static inspector method
//BpCollectionData BpCollectionInspector::inspect_bp_collection(const BpCollection&
//                                                              bp_collection) {
//
//    // data structure
//    BpCollectionData bp_coll_data;
//
//    bp_coll_data._energy_per_step = step_energies(bp_collection);
//    bp_coll_data._bend_angle_per_step = step_bend_angles(bp_collection);
//    bp_coll_data._twist_sc_per_step = step_supercoiling_twist(bp_collection);
////    bp_coll_data._force_per_bp = bp_forces(bp_collection);
////    bp_coll_data._moment_per_bp = bp_moments(bp_collection);
//
//    return bp_coll_data;
//
//};
//
//
//// step energies inspecting method
//std::vector<Real> BpCollectionInspector::step_energies(const BpCollection&
//                                                       bp_collection) {
//
//    // energy per step
//    std::vector<Real> step_energies;
//    step_energies.reserve(bp_collection.n_of_bp_steps());
//    for (Size i=0; i<bp_collection.n_of_bp_steps(); ++i) {
//        step_energies.push_back(BpCollectionElasticEnergy::
//                                single_step_elastic_energy(i, bp_collection));
//    };
//
//    return step_energies;
//
//};
//
//
//// step bending angle inspecting method
//std::vector<Real> BpCollectionInspector::step_bend_angles(const BpCollection&
//                                                          bp_collection) {
//
//    // bend angle
//    std::vector<Real> step_bend_angles;
//    step_bend_angles.reserve(bp_collection.n_of_bp_steps());
//    for (Size i=0; i<bp_collection.n_of_bp_steps(); ++i) {
//        const BpStepParams& prms(bp_collection.bp_step_params(i));
//        Real angle(prms.value(TILT)*prms.value(TILT)
//                   +prms.value(ROLL)*prms.value(ROLL));
//        step_bend_angles.push_back(std::sqrt(angle));
//    };
//
//    return step_bend_angles;
//    
//};
//
//
//// step twist of supercoiling inspecting method
//std::vector<Real>
//BpCollectionInspector::step_supercoiling_twist(const BpCollection&
//                                               bp_collection) {
//
//    // computation
//    std::vector<Vector3> tw_data =
//    BpTwistDensity::compute_twist_density(bp_collection.base_pairs());
//
//    // data filtering
//    std::vector<Real> twist_sc(tw_data.size(), FLOAT_INIT);
//    for (Size i=0; i<tw_data.size(); ++i)
//        twist_sc[i] = tw_data[i][Z];
//
//    return twist_sc;
//
//};
//
//
//// bp acting forces inspecting method
//std::vector<Vector3> BpCollectionInspector::bp_forces(const BpCollection&
//                                                      bp_collection) {
//
//    // container
//    std::vector<Vector3> forces;
//    forces.reserve(bp_collection.n_of_base_pairs());
//
//    // loop over all base pairs
//    for (Size i=0; i<bp_collection.n_of_base_pairs(); ++i) {
//
//        // first bp
//        if (i == 0) {
//            Matrix3 dmat = bp_collection.base_pair(0).column_axes();
//            dmat.transpose();
//            forces.push_back(dmat*compute_step_force(0, bp_collection));
//        }
//
//        // last bp
//        else if (i == bp_collection.n_of_base_pairs()-1) {
//            Matrix3 dmat =
//            bp_collection.base_pair(bp_collection.n_of_base_pairs()-1).
//            column_axes();
//            dmat.transpose();
//            forces.push_back(dmat*(Real(-1)*compute_step_force(0,
//                                                               bp_collection)));
//        }
//
//        // other bps
//        else {
//            Matrix3 dmat = bp_collection.base_pair(i).column_axes();
//            dmat.transpose();
//            Vector3 f1 = Real(-1)*compute_step_force(i-1, bp_collection);
//            Vector3 f2 = compute_step_force(i, bp_collection);
//            forces.push_back(dmat*(f1+f2));
//        };
//
//    };
//
//    return forces;
//
//};
//
//
//// bp acting moments inspecting method
//std::vector<Vector3> BpCollectionInspector::bp_moments(const BpCollection&
//                                                       bp_collection) {
//
//    // container
//    std::vector<Vector3> moments;
//    moments.reserve(bp_collection.n_of_base_pairs());
//
//    // loop over all base pairs
//    for (Size i=0; i<bp_collection.n_of_base_pairs(); ++i) {
//
//        // first bp
//        if (i == 0) {
//            Matrix3 dmat = bp_collection.base_pair(0).column_axes();
//            dmat.transpose();
//            moments.push_back(compute_step_moment(0, bp_collection));
//        }
//
//        // last bp
//        else if (i == bp_collection.n_of_base_pairs()-1) {
//            Matrix3 dmat =
//            bp_collection.base_pair(bp_collection.n_of_base_pairs()-1).
//            column_axes();
//            dmat.transpose();
//            moments.push_back((Real(-1)*
//                                    compute_step_moment(0, bp_collection)));
//        }
//
//        // other bps
//        else {
//            Matrix3 dmat = bp_collection.base_pair(i).column_axes();
//            dmat.transpose();
//            Vector3 m1 = Real(-1)*compute_step_moment(i-1, bp_collection);
//            Vector3 m2 = compute_step_moment(i, bp_collection);
//            moments.push_back((m1+m2));
//        };
//        
//    };
//    
//    return moments;
//    
//};
//
//
//// force computation method
//Vector3 BpCollectionInspector::compute_step_force(Size step_index,
//                                                  const BpCollection&
//                                                  bp_collection) {
//
//    // step parameters gradient
//    const BpStepParams& p = bp_collection.bp_step_params(step_index);
//    const BpStepParams& p0 =
//    bp_collection.bp_step_intrinsic_parameters(step_index);
//    const MatrixN& fmat = bp_collection.bp_step_force_constants(step_index);
//    VectorN g = BpCollectionGradient::step_parameters_gradient(p, p0, fmat);
//    Vector3 rho_g(g[3], g[4], g[5]);
//
//    // jacobian matrix T
//    const Matrix3 Tmat =
//    BpCollectionGeometry::jacobian_matrix_T(step_index, bp_collection);
//
//    // transformation
//    // this result is in kbT.A^-1
//    Vector3 force = Tmat*rho_g;
//
//    // convert to pN
//#define kBT_SCALE (4.11*10)
//    return force*kBT_SCALE;
//#undef kBT_SCALE
//
//};
//
//
//// moment computation method
//Vector3 BpCollectionInspector::compute_step_moment(Size step_index,
//                                                   const BpCollection&
//                                                   bp_collection) {
//
//    // step parameters gradient
//    const BpStepParams& p = bp_collection.bp_step_params(step_index);
//    const BpStepParams& p0 =
//    bp_collection.bp_step_intrinsic_parameters(step_index);
//    const MatrixN& fmat = bp_collection.bp_step_force_constants(step_index);
//    VectorN g = BpCollectionGradient::step_parameters_gradient(p, p0, fmat);
//    Vector3 theta_g(g[0], g[1], g[2]);
//    Vector3 rho_g(g[3], g[4], g[5]);
//
//    // sigma matrix
//    Matrix3 Sigma =
//    BpGeometryFunctions::
//    step_Sigma_matrix(BpCollectionGeometry::
//                      ZYZEulerAngles_from_step_parameters(p));
//    Sigma.transpose();
//
//    // jacobian R
//    Matrix3 Rmat = BpCollectionGeometry::jacobian_matrix_R(step_index,
//                                                           bp_collection);
//
//#define kBT_SCALE (4.11*10)
//    return kBT_SCALE*Sigma*(DEG_2_RAD*theta_g+Rmat*rho_g);
//#undef kBT_SCALE
//
//};
//
