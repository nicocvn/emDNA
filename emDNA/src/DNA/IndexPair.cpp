// IndexPair class
// Nicolas Clauvelin


#include <IndexPair.h>


// class default constructor
IndexPair::IndexPair() : m_pair_of_indexes(SizePair(Size(0), Size(0))) {};


// class constructor by copy
IndexPair::IndexPair(const IndexPair& idx_pair) :
m_pair_of_indexes(idx_pair.m_pair_of_indexes) {};


// class constructor with index initialization
IndexPair::IndexPair(Size i, Size j) : m_pair_of_indexes(SizePair(i,j)) {};


// class destructor
IndexPair::~IndexPair() {};


// copy operator
IndexPair& IndexPair::operator=(const IndexPair& idx_pair) {
    m_pair_of_indexes = idx_pair.m_pair_of_indexes;
    return *this;
};


// collection index
const Size& IndexPair::collection_index() const {
    return m_pair_of_indexes.first;
};


// local index
const Size& IndexPair::local_index() const {
    return m_pair_of_indexes.second;
};

