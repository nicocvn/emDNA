// Minimizer_Alglib abstract class


#ifndef emDNA_Minimizer_Alglib_h
#define emDNA_Minimizer_Alglib_h


#include <Alglib_Includes.h>


// enum for return code
enum AlglibMinReturnCode {
    EPSX = 0,       // relative step size is less than threshold on x
    EPSF = 1,       // relative function improvement is less than threshold on f
    EPSG = 2,       // relative gradient norm is less than threshold on g
    MAXIT = 3,      // maximum number of iterations is reached
    BADCONDS = 4,   // conditions are too stringent
    BADGRAD = 5,    // error in gradient computation (for debug purposes only)
    UNKNOWN = 6     // unknown return code
};
inline std::ostream & operator<<(std::ostream& os, AlglibMinReturnCode R) {
    switch (R) {
        case EPSX:
            return (os << "EPSX");
        case EPSF:
            return (os << "EPSF");
        case EPSG:
            return (os << "EPSG");
        case MAXIT:
            return (os << "MAXIT");
        case BADCONDS:
            return (os << "BADCONDS");
        case BADGRAD:
            return (os << "BADGRAD");
        case UNKNOWN:
            return (os << "UNKNOWN");
        default:
            return (os << ((Integer)R));
    };
};


// alglib minimizer results data structure
struct AlglibMinResults {
    AlglibMinReturnCode _return_code;
    Size _iterations;
    Size _f_evaluations;
    Real _initial_function;
    Real _optimal_function;
    VectorN _optimal_point;
    VectorN _optimal_gradient;
};


// pointer typedef for the objective function
typedef void (*ObjFunction)(const alglib::real_1d_array &x,
                            double &func,
                            alglib::real_1d_array &grad,
                            void *ptr);
typedef void (*CallbackFunction)(const alglib::real_1d_array& x,
                                 double func,
                                 void *data_ptr);


// Minimizer_Alglib class
class Minimizer_Alglib {


public:

    // virtual destructor
    virtual ~Minimizer_Alglib() {};

    // minimizer function modifier
    void set_minimizer_function(ObjFunction f) {
        m_objective_function = f;
    };

    // minimizer starting point modifier
    virtual void set_minimizer_starting_point(const VectorN& start_x) {
        m_minimizer_x = real_1d_array_with_size(start_x.size());
        VectorN_TO_real_1d_array(start_x, m_minimizer_x);
    };

    // minimizer settings accessor/modifier
    virtual AlglibMinSettings minimizer_settings() const =0;
    virtual void set_minimizer_settings(const AlglibMinSettings&
                                        minim_settings) =0;

    // minimizer scaling modifier
    virtual void set_minimizer_x_scalings(const VectorN& scale_x) =0;

    // minimizer boundaries modifier
    // should only be implemented for BLEIC algorithm
    virtual void set_minimizer_boundaries(const VectorN& lower_bc,
                                          const VectorN& upper_bc) =0;

    // minimize methods
    virtual AlglibMinResults minimize(CallbackFunction callback_f,
                                      void* data_ptr) =0;

    // minimizer gradient checking method
    virtual AlglibMinResults perform_gradient_checking(Real step_size,
                                                       Integer* idx,
                                                       void* data_ptr) =0;


protected:

    // private constructors and copy operator
    Minimizer_Alglib() :
        m_objective_function(nullptr), m_minimizer_x(), m_scales() {};
    Minimizer_Alglib(const Minimizer_Alglib& alglib_minimizer) :
        m_objective_function(alglib_minimizer.m_objective_function),
        m_minimizer_x(alglib_minimizer.m_minimizer_x),
        m_scales(alglib_minimizer.m_scales) {};
    Minimizer_Alglib& operator=(const Minimizer_Alglib& alglib_minimizer) {
        m_objective_function = alglib_minimizer.m_objective_function;
        m_minimizer_x = alglib_minimizer.m_minimizer_x;
        m_scales = alglib_minimizer.m_scales;
        return *this;
    };

    // alglib return code translation
    template <class MinimizerType>
    static AlglibMinReturnCode
    return_code_translation(const MinimizerType& rep) {
        alglib::ae_int_t ret = rep.terminationtype;
        if (ret == -7)
            return BADGRAD;
        if (ret == -3)
            return BADCONDS;
        else if (ret == 1)
            return EPSF;
        else if (ret == 2)
            return EPSX;
        else if (ret == 4)
            return EPSG;
        else if (ret == 5)
            return MAXIT;
        else if (ret == 7)
            return BADCONDS;
        else
            return UNKNOWN;
    };

    // alglib minimizer data
    ObjFunction m_objective_function;
    alglib::real_1d_array m_minimizer_x;
    VectorN m_scales;
    
    
};


#endif  // emDNA_Minimizer_Alglib_h
