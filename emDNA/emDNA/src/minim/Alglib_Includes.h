// Alglib Includes
// Nicolas Clauvelin


#ifndef emDNA_Alglib_Includes_h
#define emDNA_Alglib_Includes_h


#include <emDNA_Includes.h>


// alglib includes
#include <stdafx.h>
#include <optimization.h>


// vector conversion functions
alglib::real_1d_array real_1d_array_with_size(Size n);
void real_1d_array_TO_VectorN(const alglib::real_1d_array& source,
                              VectorN& target);
void VectorN_TO_real_1d_array(const VectorN& source,
                              alglib::real_1d_array& target);


#endif  // emDNA_Alglib_Includes_h
