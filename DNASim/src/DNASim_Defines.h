// DNASim Defines header
// Nicolas Clauvelin


// project-wide header file containing typedefs, defines and useful functions


// DNASim precision
#define DNASIM_DOUBLE_PREC


// DNASim numerical types
using Size = std::size_t;
using UInteger = unsigned int;
using Integer = int;
#ifdef DNASIM_SINGLE_PREC
    using Real = float;
#endif	// DNASIM_SINGLE_PREC
#ifdef DNASIM_DOUBLE_PREC
    using Real = double;
#endif	// DNASIM_DOUBLE_PREC


// numerical constants
#ifdef DNASIM_SINGLE_PREC
    #define ZERO_EPS Real(1E-15)        // the machine epsilon is (guess) 1e-16
	#define ZERO_TOL Real(1E-6)
	#define ZERO_TOL_SQ Real(1E-12)
#endif	// DNASIM_SINGLE_PREC
#ifdef DNASIM_DOUBLE_PREC
    #define ZERO_EPS Real(1E-15)        // the machine epsilon is (guess) 1e-16
    #define ZERO_TOL Real(1E-10)
    #define ZERO_TOL_SQ Real(1E-20)
#endif	// DNASIM_DOUBLE_PREC
#define FLOAT_INIT Real(0.0)
#define F_PI Real(M_PI)
#define DEG_2_RAD F_PI/Real(180.0)
#define RAD_2_DEG Real(180.0)/F_PI


// real number width for output
#ifdef DNASIM_SINGLE_PREC
	#define REAL_WIDTH 6
#endif
#ifdef DNASIM_DOUBLE_PREC
	#define REAL_WIDTH 6 // changed from 10 to 6 (11 Dec 2018)
#endif


namespace DNASim {


    // types renaming
    using KeyValueData = std::map<std::string,std::string>;


    // enum class for PRN random seeding
    enum class RandomSeed : bool {
        Yes = true, No = false
    };


    // pair of real numbers for complex number storage
    using RealPair = std::pair<Real,Real>;

        
    // array structure
    template <class ArrayType> struct Array {

        // array size
        Size _size;

        // array values
        const ArrayType* const _values;

        // constructor
        Array(Size n, const ArrayType* const values) :
            _size(n), _values(values) {};

    };


    // STL vector output template function
    template <class T>
    void print_vector(const std::vector<T>& data,
                      const char& opening = '{',
                      const char& closing = '}',
                      const char& separator = ',') {

        Size n = data.size();
        if (n == 0)
            return;

        std::cout << opening;
        for (Size i=0; i<n; ++i) {
            std::cout << data[i];
            if (i != n-1)
                std::cout << separator;
        };
        std::cout << closing;
        
    };


    // pointer delete function
	template <class PointerType> void pointer_delete(PointerType* p) {
		delete p;		// implementation-safe
		p = nullptr;
	};


}

// DNASim SVN revision
#define DNASIM_CODE_REV "1114:1127"
