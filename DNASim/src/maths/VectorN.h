// VectorN class header file
// Nicolas Clauvelin

// n-dimensional Real vector implementation based on Eigen


#ifndef DNASim_VectorN_h
#define DNASim_VectorN_h


#include <Eigen_Includes.h>


namespace DNASim {


    class VectorN {


    public:

        // constructors
		// default constructor creates an empty vector
		VectorN();
        VectorN(const VectorN& v) = default;
        VectorN(VectorN&& v);
		VectorN(Size n, Real val);
		VectorN(const Array<Real>& array);
        VectorN(const std::initializer_list<Real>& v);
        explicit VectorN(const std::string& s);
        explicit VectorN(const std::vector<Real>& stl_vector);
        ~VectorN() = default;

        // size accessor/modifier
        // all entries are set to zero upon resizing
		Size size() const;
		void resize(Size n);

        // component const-accessor
        inline const Real&	operator[](Size i) const {
            return m_vec[static_cast<Integer>(i)];
        };

        // component accessor
        inline Real& operator[](Size i) {
            return m_vec[static_cast<Integer>(i)];
        };

		// vector manipulations
        void clear();
		Real norm() const;
		void normalize();
        void diagonal_scale(const VectorN& v);
        Real dot(const VectorN& v) const;

        // vector slicing
        // slice() returns the vector [first,last)
        VectorN slice(Size first, Size last) const;
        void set_slice(Size first, const VectorN& slice);

        // pack values methods
        Array<Real> pack_values() const;
        std::vector<Real> std_vector() const;

		// operators
        VectorN& operator=(const VectorN& v) = default;
        VectorN& operator=(VectorN&& v);
		VectorN& operator+=(const VectorN& v);
		VectorN& operator-=(const VectorN& v);
		VectorN& operator*=(Real a);

        // vector operators
		friend VectorN operator+(const VectorN& v1, const VectorN& v2);
        friend VectorN operator-(const VectorN& v1, const VectorN& v2);
        friend VectorN operator*(const VectorN& v, const Real& a);
        friend VectorN operator*(const Real& a, const VectorN& v);

        // output operator
        friend std::ostream& operator<<(std::ostream& os, const VectorN& v);


    private:

        // private constructor with initialization from an Eigen object
        explicit VectorN(const EigenVectorX vec);

        EigenVectorX m_vec;


    };


}



#endif  // DNASim_VectorN_h
