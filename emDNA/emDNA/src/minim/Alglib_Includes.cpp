// Alglib Includes
// Nicolas Clauvelin


#include "Alglib_Includes.h"


// sized vector function
alglib::real_1d_array real_1d_array_with_size(Size n) {
    alglib::real_1d_array v;
    v.setlength((alglib::ae_int_t)n);
    return v;
};


// real_1d_array -> VectorN
void real_1d_array_TO_VectorN(const alglib::real_1d_array& source,
                              VectorN& target) {

    // size
    Size n_elems = (Size)source.length();

    // filling
    for (Size i=0; i<n_elems; ++i)
        target[i] = source[(alglib::ae_int_t)i];

};


// VectorN -> real_1d_array
void VectorN_TO_real_1d_array(const VectorN& source,
                              alglib::real_1d_array& target) {

    // size
    alglib::ae_int_t n_elems = (alglib::ae_int_t)source.size();

    // filling
    for (alglib::ae_int_t i=0; i<n_elems; ++i)
        target[i] = source[(Size)i];

};

