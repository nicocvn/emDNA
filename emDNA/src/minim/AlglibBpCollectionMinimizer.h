// AlglibBpCollectionMinimizer class
// Nicolas Clauvelin


// alglib gradient-based minimization class


#ifndef emDNA_AlglibBpCollectionMinimizer_h
#define emDNA_AlglibBpCollectionMinimizer_h


#include <Minimizer_Alglib.h>


// minimization algorithm enum type
enum MinimizationAlgo {
    LBFGS_ALGO = 0,
    CG_ALGO = 1
};

// data structure for callback functions
struct CallbackData {
    std::shared_ptr<BpCollection_Interface> _bp_coll_ptr;
    CallbackOptions _opts;
    Size _callback_counter;
};


class AlglibBpCollectionMinimizer {


public:

    // static alglib minimization method
    static AlglibMinResults
    optimize_bp_collection(const std::shared_ptr<BpCollection_Interface>&
                           interface_ptr,
                           const AlglibMinSettings& minimization_settings,
                           const MinimizationAlgo& minim_algo,
                           const CallbackOptions& callback_opts);


private:

    // minimizer instantiation method from algorithm type
    static
    std::unique_ptr<Minimizer_Alglib> create_minimizer(const MinimizationAlgo&
                                                       minim_algo);

    // private constructors
    AlglibBpCollectionMinimizer() = delete;
    ~AlglibBpCollectionMinimizer() = delete;

    // deleted copy/move constructors and operators
    AlglibBpCollectionMinimizer(const AlglibBpCollectionMinimizer&
                                alglib_minimizer) = delete;
    AlglibBpCollectionMinimizer(AlglibBpCollectionMinimizer&&
                                alglib_minimizer) = delete;
    AlglibBpCollectionMinimizer& operator=(const AlglibBpCollectionMinimizer&
                                           alglib_minimizer) = delete;
    AlglibBpCollectionMinimizer& operator=(AlglibBpCollectionMinimizer&&
                                           alglib_minimizer) = delete;


};


// private callback function
void callback_function(const alglib::real_1d_array& x,
                       double func,
                       void *data_ptr);


// alglib objective function wrapper
void alglib_obj_function(const alglib::real_1d_array& x,
                         double& func, alglib::real_1d_array& grad,
                         void *data_ptr);


#endif  // emDNA_AlglibBpCollectionMinimizer_h
