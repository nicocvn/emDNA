// AlglibBpCollectionMinimizer class
// Nicolas Clauvelin


#include "DNA/BpCollection_Interface.h"
#include "minim/Minimizer_AlglibCG.h"
#include "minim/Minimizer_AlglibLBFGS.h"
#include "minim/AlglibBpCollectionMinimizer.h"


// static alglib minimization method
AlglibMinResults AlglibBpCollectionMinimizer::
optimize_bp_collection(const std::shared_ptr<BpCollection_Interface>&
                       interface_ptr,
                       const AlglibMinSettings& minimization_settings,
                       const MinimizationAlgo& minim_algo,
                       const CallbackOptions& callback_opts) {

    // initial energy
    const Real E0 = interface_ptr->bp_collection_energy();

    // new alglib minimizer
    std::unique_ptr<Minimizer_Alglib> minimizer = create_minimizer(minim_algo);

    // starting point
    minimizer->set_minimizer_starting_point(interface_ptr->
                                            bp_collection_inline_free_dofs());

    // objective function
    minimizer->set_minimizer_function(alglib_obj_function);

    // settings
    minimizer->set_minimizer_settings(minimization_settings);

    // scalings
    minimizer->set_minimizer_x_scalings(interface_ptr->
                                        bp_collection_free_dofs_scalings());

    // callback data structure
    CallbackData callback_data;
    callback_data._bp_coll_ptr = interface_ptr;
    callback_data._opts = callback_opts;
    callback_data._callback_counter = 0;

    // clean obsolete file
    if (callback_opts._output_collection_frequency != 0)
        remove("tmp_confs.txt");

    // optimization
    AlglibMinResults minim_results = minimizer->minimize(&callback_function,
                                                         (void*)&callback_data);

    // clean progress indicate
    if (callback_opts._output_current_energy) {
        std::cout << "\r\n\n";
        std::cout.flush();
    };

    // add the initial energy to the results
    minim_results._initial_function = E0;
    
    return minim_results;

};


// minimizer instantiation method from algorithm type
std::unique_ptr<Minimizer_Alglib>
AlglibBpCollectionMinimizer::create_minimizer(const MinimizationAlgo&
                                              minim_algo) {

    std::unique_ptr<Minimizer_Alglib> minimizer_ptr = nullptr;

    switch (minim_algo) {
        case LBFGS_ALGO:
            minimizer_ptr =
            std::unique_ptr<Minimizer_AlglibLBFGS>(new Minimizer_AlglibLBFGS());
            break;
        case CG_ALGO:
            minimizer_ptr =
            std::unique_ptr<Minimizer_AlglibCG>(new Minimizer_AlglibCG());
            break;
        default:
            minimizer_ptr = nullptr;
            DS_ASSERT(minimizer_ptr != nullptr,
                      "Alglib minimizer initialization with a non-supported "
                      "minimization algorithm");
            break;
    };

    return minimizer_ptr;

};


// private callback function
void callback_function(const alglib::real_1d_array& x,
                       double func,
                       void *data_ptr) {

    // static cast
    CallbackData* callback_data = static_cast<CallbackData*>(data_ptr);

    // exit if no callback output
    if(!callback_data->_opts._callback_output)
        return;

    // print current energy value
    if (callback_data->_opts._output_current_energy) {
        std::cout << std::setprecision(REAL_WIDTH);
        std::cout << "  ### current energy = " << func << "\r";
        std::cout.flush();
    };

    // configuration output
    if (callback_data->_opts._output_collection_frequency != 0) {

        // counter increment
        ++callback_data->_callback_counter;

        // limit reached
        if (callback_data->_callback_counter ==
            callback_data->_opts._output_collection_frequency) {

            // bp collection
            const BpCollection bp_coll =
            callback_data->_bp_coll_ptr->bp_collection();

            // bp collection inline formatting
            std::stringstream bp_stream;
            const Size n_bps = bp_coll.n_of_base_pairs();
            bp_stream << "{";
            for (Size i=0; i<n_bps-1; ++i)
                bp_stream << bp_coll.base_pair(i) << ",";
            bp_stream << bp_coll.base_pair(n_bps-1) << "}";

            // output
            OutputFileHandler output_file("tmp_confs.txt", ".");
            output_file.open_in_append_mode();
            output_file.write_line(bp_stream.str());
            output_file.close();

            // reset counter
            callback_data->_callback_counter = 0;

        };

    };

};


// alglib objective function wrapper
void alglib_obj_function(const alglib::real_1d_array& x,
                         double& func, alglib::real_1d_array& grad,
                         void *data_ptr) {

    // pointer checking
    DS_ASSERT(data_ptr != NULL,
              "uninitialized data pointer in alglib objective function");

    // static cast
    CallbackData *callback_data = static_cast<CallbackData*>(data_ptr);

    // bp collection interface
    std::shared_ptr<BpCollection_Interface> bp_collection_interface =
    callback_data->_bp_coll_ptr;

    // convert dofs from alglib vector to VectorN
    VectorN inline_free_dofs((Size)x.length(), FLOAT_INIT);
    real_1d_array_TO_VectorN(x, inline_free_dofs);

    // bp collection update with inline dofs
    // these dofs correspond to the free dofs
    bp_collection_interface->update_with_new_inline_free_dofs(inline_free_dofs);

    // energy
    func = bp_collection_interface->bp_collection_energy();

    // gradient
    VectorN gradient(bp_collection_interface->
                     bp_collection_free_dofs_gradient());
    VectorN_TO_real_1d_array(gradient, grad);

    // check gradient dimensions
    DS_ASSERT(gradient.size() == inline_free_dofs.size(),
              "wrong gradient size compared to free dofs size in "
              "alglig objective function");
    DS_ASSERT((Size)grad.length() == inline_free_dofs.size(),
              "wrong alglib gradient size compared to free dofs size in "
              "alglig objective function");

};


