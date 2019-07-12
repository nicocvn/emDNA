// VectorN class implementation file
// Nicolas Clauvelin


#include <EnhancedString.h>
#include <VectorN.h>


namespace DNASim {


    // default constructor
    // initializes the vector as empty
    VectorN::VectorN() : m_vec() {};


    // constructor by moving
    VectorN::VectorN(VectorN&& v) : VectorN() {
        // Eigen types are not movable at this point so we rely on the swap
        m_vec.swap(v.m_vec);
    };


    // constructor with size and default value initialization
    // all the components of the vector are initialized to val
    VectorN::VectorN(Size n, Real val) :
    m_vec(EigenVectorX::Constant((Integer)n, val)) {};


    // constructor with initialization from an Array
    VectorN::VectorN(const Array<Real>& array) :
    m_vec(Eigen::Map<const EigenVectorX>(&(array._values[0]),
                                         (Integer)array._size)) {};


    // constructor from a initializer list
    VectorN::VectorN(const std::initializer_list<Real>& v) :
    VectorN(v.size(), FLOAT_INIT) {
        for (Size i=0, end=v.size(); i<end; ++i)
            m_vec[i] = *(v.begin()+i);
    };


    // constructor with initialization from a string
    // the string is expected to be of the form "{x1, x2, x3, ..., xN}" where N
    // is the size of the vector
    VectorN::VectorN(const std::string& s) : VectorN() {

        // string cleaning
		EnhancedString es(s);
		es.erase_blank_characters();
		es.erase_characters({'{','}'});

		// tokenizing
		std::vector<EnhancedString> tokens;
		es.tokenize(tokens, ',');

		// vector size
		Size n = tokens.size();
		resize(n);

		// values assignement
		for (Size i=0; i<n; i++)
			m_vec[(Integer)i] = tokens[i].convert<Real>();

    };


    // constructor with initialization from a std::vector
    VectorN::VectorN(const std::vector<Real>& stl_vector) :
    m_vec(Eigen::Map<const EigenVectorX>(&stl_vector[0],
                                         (Integer)stl_vector.size())) {};


    // size accessor
    Size VectorN::size() const { return (Size)m_vec.size(); };

    // size modifier
    // all the components are set to zero upon resizing
    void VectorN::resize(Size n) {
        m_vec.resize((Integer)n);
        clear();
    };


    // clear method
    // set all the components to zero
    void VectorN::clear() { m_vec.setZero(); };

    /// norm method
    // return the Euclidean norm of the vector
    Real VectorN::norm() const { return m_vec.norm(); };

    // normalize method
    void VectorN::normalize() { m_vec.normalize(); };

    // diagonal scale method
    void VectorN::diagonal_scale(const VectorN& v) {
        EigenMatrixX m(m_vec.asDiagonal());
        m_vec = m*v.m_vec;
    };

    // scalar product computation method
    Real VectorN::dot(const VectorN& v) const {
        return m_vec.dot(v.m_vec);
    };


    // vector slicing methods
    // returns the vector [first,last)
    VectorN VectorN::slice(Size first, Size last) const {
        return VectorN(m_vec.segment(first, last-first));
    };
    void VectorN::set_slice(Size first, const VectorN& slice) {
        m_vec.segment(first, slice.size()) = slice.m_vec;
    };


    // pack values
    Array<Real> VectorN::pack_values() const {
        return Array<Real>((Size)m_vec.size(), &(m_vec(0)));
    };
    std::vector<Real> VectorN::std_vector() const {
        return std::vector<Real>(m_vec.data(), m_vec.data()+m_vec.size());
    };


    // operators
    VectorN& VectorN::operator=(VectorN&& v) {
        // Eigen types are not movable at this point so we rely on the swap
        m_vec.swap(v.m_vec);
        return *this;
    };
    VectorN& VectorN::operator+=(const VectorN& v) {
        m_vec += v.m_vec;
        return *this;
    };
    VectorN& VectorN::operator-=(const VectorN& v)  {
        m_vec -= v.m_vec;
        return *this;
    };
    VectorN& VectorN::operator*=(Real a) {
        m_vec *= a;
        return *this;
    };

    
    // vector operators
    VectorN operator+(const VectorN& v1, const VectorN& v2) {
        return VectorN(v1.m_vec+v2.m_vec);
    };
    VectorN operator-(const VectorN& v1, const VectorN& v2) {
        return VectorN(v1.m_vec-v2.m_vec);
    };
    VectorN operator*(const VectorN& v, const Real& a) {
        return VectorN(a*v.m_vec);
    };
    VectorN operator*(const Real& a, const VectorN& v) {
        return VectorN(a*v.m_vec);
    };


    // output operator
    std::ostream& operator<<(std::ostream& os, const VectorN& v) {
        os << '{' << std::fixed << std::setprecision(REAL_WIDTH);
        for (Size i=0, end=v.size()-1; i<end; ++i)
            os << v[i] << ", ";
        os << v[v.size()-1];
		os << '}';
		return os;
    };


    // private constructor with initialization from an Eigen object
    VectorN::VectorN(const EigenVectorX vec) : m_vec(vec) {};
    

}
