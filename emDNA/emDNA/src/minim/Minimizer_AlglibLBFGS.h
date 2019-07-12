// Minimizer_AlglibLBFGS class
// Nicolas Clauvelin


// conjugate gradient minimization implementation based on alglib::minlbfgs
// see  http://www.alglib.net/translator/man/manual.cpp.html#unit_minlbfgs

// ** stopping conditions (see AlglibLBGFSSettings):
//
// step size threshold:
// optimization finishes if there is the step size is less than the threshold
// |dx| < _threshold_x
//
// function improvement:
// optimization finishes if the following condition is satisfied at the (k+1)-th
// iteration:
// |f(k+1)-f(k)|<=_threshold_f*max{|f(k)|,|f(k+1)|,1}
//
// gradient norm:
// optimization finishes if the norm of the gradient is less than the threshold:
// |g| < _threshold_g
//
// if all thresholds are set to zero the stopping criterion is based on relative
// step size


// ** usage:
//
// 1- set the starting point -> set_minimizer_starting_point(...)
// 2- set the objective function -> set_minimizer_function(...)
// 3- set the minimizer settings -> set_minimizer_settings(...)
//    eventually set the variable scales -> set_minimizer_x_scalings(...)
// 4- perform minimization -> minimize()
// 5- if needed restart (settings can be changed between restarts)


// ** objective function format:
//
// func(const alglib::real_1d_array &x,         // point
//      double &func,                           // objective value
//      alglib::real_1d_array &grad,            // objective gradient
//      void *ptr)                              // parameters pointer


#ifndef emDNA_Minimizer_AlglibLBFGS_h
#define emDNA_Minimizer_AlglibLBFGS_h


#include <emDNA_Includes.h>
#include <Minimizer_Alglib.h>


// Minimizer_AlglibLBFGS class
class Minimizer_AlglibLBFGS : public Minimizer_Alglib {


public:

    // constructors
    Minimizer_AlglibLBFGS();
    Minimizer_AlglibLBFGS(const Minimizer_AlglibLBFGS& alglib_minimizer);
    ~Minimizer_AlglibLBFGS();

    // copy operator
    Minimizer_AlglibLBFGS& operator=(const Minimizer_AlglibLBFGS&
                                     alglib_minimizer);


    // minimizer starting point modifier
    void set_minimizer_starting_point(const VectorN& start_x);

    // minimizer settings accessor/modifier
    AlglibMinSettings minimizer_settings() const;
    void set_minimizer_settings(const AlglibMinSettings& minim_settings);

    // minimzer scaling modifier
    void set_minimizer_x_scalings(const VectorN& scale_x);

    // minimizer boundaries modifier
    // should only be implemented for BLEIC algorithm
    void set_minimizer_boundaries(const VectorN& lower_bc,
                                  const VectorN& upper_bc) {
        DS_ASSERT(0==1,
                  "setting boundaries in non-BLEIC minimizer");
    };

    // minimize methods
    AlglibMinResults minimize(CallbackFunction callback_f,
                              void* data_ptr);

    // minimizer gradient checking method
    AlglibMinResults perform_gradient_checking(Real step_size,
                                               Integer* idx,
                                               void* data_ptr);


private:

    // Alglib minimizer data
    alglib::minlbfgsstate m_minimizer_state;
    alglib::minlbfgsreport m_minimizer_report;


};


#endif  // emDNA_Minimizer_AlglibLBFGS_h
