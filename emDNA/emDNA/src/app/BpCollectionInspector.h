//// BpCollectionInspector class
//// Nicolas Clauvelin
//
//
//#ifndef emDNA_BpCollectionInspector_h
//#define emDNA_BpCollectionInspector_h
//
//
//#include <emDNA_Includes.h>
//class BpCollection;
//
//
//// data structure for bp collection data
//struct BpCollectionData {
//    std::vector<Real> _energy_per_step;
//    std::vector<Real> _bend_angle_per_step;
//    std::vector<Real> _twist_sc_per_step;
//};
//
//
//class BpCollectionInspector {
//
//
//public:
//
//    // static output to file method
//    static void write_inspector_data(const std::string& output_file,
//                                     const BpCollection& bp_collection);
//
//    // static inspector method
//    static BpCollectionData inspect_bp_collection(const BpCollection&
//                                                  bp_collection);
//
//
//private:
//
//    // private inspector methods
//    static std::vector<Real> step_energies(const BpCollection& bp_collection);
//    static std::vector<Real> step_bend_angles(const BpCollection&
//                                              bp_collection);
//    static std::vector<Real> step_supercoiling_twist(const BpCollection&
//                                                     bp_collection);
//    static std::vector<Vector3> bp_forces(const BpCollection& bp_collection);
//    static std::vector<Vector3> bp_moments(const BpCollection& bp_collection);
//
//    // force computation method
//    static Vector3 compute_step_force(Size step_index,
//                                      const BpCollection& bp_collection);
//
//    // moment computation method
//    static Vector3 compute_step_moment(Size step_index,
//                                       const BpCollection& bp_collection);
//
//    // private constructors and copy operator
//    BpCollectionInspector();
//    BpCollectionInspector(const BpCollectionInspector& bp_collection_inspector);
//    BpCollectionInspector& operator=(const BpCollectionInspector&
//                                     bp_collection_inspector);
//
//
//};
//
//
//#endif  // emDNA_BpCollectionInspector_h
