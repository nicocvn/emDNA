// Minimizer_AlglibLBFGS class
// Nicolas Clauvelin


#include <Minimizer_AlglibLBFGS.h>


// class default constructor
Minimizer_AlglibLBFGS::Minimizer_AlglibLBFGS() :
Minimizer_Alglib::Minimizer_Alglib(),
m_minimizer_state(),
m_minimizer_report() {};


// class constructor by copy
Minimizer_AlglibLBFGS::Minimizer_AlglibLBFGS(const Minimizer_AlglibLBFGS&
                                             alglib_minimizer) :
Minimizer_Alglib::Minimizer_Alglib(alglib_minimizer),
m_minimizer_state(alglib_minimizer.m_minimizer_state),
m_minimizer_report(alglib_minimizer.m_minimizer_report) {};


// class destructor
Minimizer_AlglibLBFGS::~Minimizer_AlglibLBFGS() {};


// copy operator
Minimizer_AlglibLBFGS&
Minimizer_AlglibLBFGS::operator=(const Minimizer_AlglibLBFGS& alglib_minimizer) {
    Minimizer_Alglib::operator=(alglib_minimizer);
    m_minimizer_state = alglib_minimizer.m_minimizer_state;
    m_minimizer_report = alglib_minimizer.m_minimizer_report;
    return *this;
};


// minimizer starting point modifier
#define LBFGS_M_PARAM 7
void Minimizer_AlglibLBFGS::set_minimizer_starting_point(const VectorN& start_x) {
    Minimizer_Alglib::set_minimizer_starting_point(start_x);
    alglib::minlbfgscreate(LBFGS_M_PARAM, m_minimizer_x, m_minimizer_state);
};
#undef LBFGS_M_PARAM


// minimizer settings accessor
AlglibMinSettings Minimizer_AlglibLBFGS::minimizer_settings() const {

    AlglibMinSettings minim_settings;

    // max iterations
    minim_settings._max_iterations = (Size)m_minimizer_state.c_ptr()->maxits;

    // step size threshold
    minim_settings._threshold_dx = m_minimizer_state.c_ptr()->epsx;

    // function threshold
    minim_settings._threshold_f = m_minimizer_state.c_ptr()->epsf;

    // gradient threshold
    minim_settings._threshold_g = m_minimizer_state.c_ptr()->epsg;

    // max step size
    minim_settings._max_step_size = m_minimizer_state.c_ptr()->stpmax;

    return minim_settings;

};


// minimizer settings modifier
void Minimizer_AlglibLBFGS::set_minimizer_settings(const AlglibMinSettings&
                                                   minim_settings) {

    // minLBFGS settings
    alglib::minlbfgssetcond(m_minimizer_state,
                            minim_settings._threshold_g,
                            minim_settings._threshold_f,
                            minim_settings._threshold_dx,
                            (alglib::ae_int_t)minim_settings._max_iterations);

    // upper bound on the step size
    if (minim_settings._max_step_size != FLOAT_INIT)
        alglib::minlbfgssetstpmax(m_minimizer_state,
                                  minim_settings._max_step_size);

};


// minimzer scaling modifier
void Minimizer_AlglibLBFGS::set_minimizer_x_scalings(const VectorN& scale_x) {
    alglib::real_1d_array scalings(real_1d_array_with_size(scale_x.size()));
    VectorN_TO_real_1d_array(scale_x, scalings);
    alglib::minlbfgssetscale(m_minimizer_state, scalings);
    alglib::minlbfgssetprecscale(m_minimizer_state);
    m_scales = scale_x;
};


// minimize method
AlglibMinResults Minimizer_AlglibLBFGS::minimize(CallbackFunction callback_f,
                                                 void* data_ptr) {

    // turn on reporting
    // needed to use callback function
    alglib::minlbfgssetxrep(m_minimizer_state, true);

    // optimization
    // it is possible to pass a callback for each iterations (see doc)
    // prototype for callback function
    // void  (*rep)(const real_1d_array &x, double func, void *ptr)
    alglib::minlbfgsoptimize(m_minimizer_state, m_objective_function,
                             callback_f,
                             data_ptr);

    // results
    alglib::minlbfgsresults(m_minimizer_state,
                            m_minimizer_x,
                            m_minimizer_report);
    VectorN optimal_x((Size)m_minimizer_x.length(), FLOAT_INIT);
    VectorN optimal_g((Size)m_minimizer_x.length(), FLOAT_INIT);
    real_1d_array_TO_VectorN(m_minimizer_x, optimal_x);
    real_1d_array_TO_VectorN(m_minimizer_state.g, optimal_g);

    // structure data
    AlglibMinResults min_results;
    min_results._return_code =
    Minimizer_Alglib::return_code_translation(m_minimizer_report);
    min_results._iterations = (Size)m_minimizer_report.iterationscount;
    min_results._f_evaluations = (Size)m_minimizer_report.nfev;
    min_results._optimal_function = m_minimizer_state.f;
    min_results._optimal_point = optimal_x;
    min_results._optimal_gradient = optimal_g;

    return min_results;

};


// minimizer gradient checking method
AlglibMinResults
Minimizer_AlglibLBFGS::perform_gradient_checking(Real step_size,
                                                 Integer* idx,
                                                 void* data_ptr) {

    // minimization settings
    AlglibMinSettings settings = minimizer_settings();
    settings._max_iterations = 1;
    set_minimizer_settings(settings);

    // enable gradient checking
    alglib::minlbfgssetgradientcheck(m_minimizer_state, step_size);

    // optimization
    // it is possible to pass a callback for each iterations (see doc)
    alglib::minlbfgsoptimize(m_minimizer_state, m_objective_function,
                             NULL, data_ptr);

    // results
    alglib::minlbfgsresults(m_minimizer_state,
                            m_minimizer_x,
                            m_minimizer_report);
    VectorN optimal_x((Size)m_minimizer_x.length(), FLOAT_INIT);
    VectorN optimal_g((Size)m_minimizer_x.length(), FLOAT_INIT);
    real_1d_array_TO_VectorN(m_minimizer_x, optimal_x);
    real_1d_array_TO_VectorN(m_minimizer_state.g, optimal_g);

    // index
    *idx = m_minimizer_report.varidx;

    // structure data
    AlglibMinResults min_results;
    min_results._return_code =
    Minimizer_Alglib::return_code_translation(m_minimizer_report);
    min_results._iterations = (Size)m_minimizer_report.iterationscount;
    min_results._f_evaluations = (Size)m_minimizer_report.nfev;
    min_results._optimal_function = m_minimizer_state.f;
    min_results._optimal_point = optimal_x;
    min_results._optimal_gradient = optimal_g;

    return min_results;

};
