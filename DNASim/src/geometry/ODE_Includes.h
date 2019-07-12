// ODE Includes
// Nicolas Clauvelin


#ifndef DNASim_ODE_Includes_h
#define DNASim_ODE_Includes_h


#include <DNASim_Includes.h>


#ifdef WITH_ODE_COLLISION

// suppress warnings
#pragma clang diagnostic ignored "-Weffc++"
#pragma clang diagnostic ignored "-Wwrite-strings"
#pragma clang diagnostic ignored "-Wdeprecated-writable-strings"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wwrite-strings"

#ifdef DNASIM_SINGLE_PREC
    #define dSINGLE
    #define dIDESINGLE
#endif	// DNASIM_SINGLE_PREC
#ifdef DNASIM_DOUBLE_PREC
    #define dDOUBLE
    #define dIDEDOUBLE
#endif	// DNASIM_DOUBLE_PREC

#include <ode/ode.h>

#endif	// WITH_ODE_COLLISION


#endif  // DNASim_ODE_Includes_h
