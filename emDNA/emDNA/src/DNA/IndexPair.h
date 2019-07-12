// IndexPair class
// Nicolas Clauvelin


// class used for symbolic numbering of bp steps within the bp collection

// the collection index is the index of the item within the bp collection
// the local index is a secondary index numbering the same item within a subset


#ifndef emDNA_IndexPair_h
#define emDNA_IndexPair_h


#include <emDNA_Includes.h>


class IndexPair {


public:

    // constructors
    IndexPair();
    IndexPair(const IndexPair& idx_pair);
    IndexPair(Size i, Size j);
    ~IndexPair();

    // copy operator
    IndexPair& operator=(const IndexPair& idx_pair);

    // collection index
    const Size& collection_index() const;

    // local index
    const Size& local_index() const;


private:

    // pair of indexes
    SizePair m_pair_of_indexes;


};


#endif  // emDNA_IndexPair_h
