// SparseIndex class
// Nicolas Clauvelin


#include <SparseIndex.h>


// constructor with range initialization
SparseIndex::SparseIndex(const Size& index_start, const Size& index_end) :
m_indices(index_end-index_start, IndexSlot(Size(0),IndexType::Clean)) {

    // index values
    Size current_val = index_start;
    for (auto it = m_indices.begin(); it != m_indices.end(); ++it) {
        it->first = current_val;
        ++current_val;
    };

};


// number of clean indices accessor
Size SparseIndex::n_clean_indices() const {
    return (Size)std::count_if(m_indices.begin(), m_indices.end(),
                               [](const IndexSlot& idx_slot) -> bool {
                                   return idx_slot.second == IndexType::Clean;
                               });
};


// number of flagged indices accessor
Size SparseIndex::n_flagged_indices() const {
    return (Size)std::count_if(m_indices.begin(), m_indices.end(),
                               [](const IndexSlot& idx_slot) -> bool {
                                   return idx_slot.second == IndexType::Flagged;
                               });
};


// ---
// const_iterator class

// constructor
SparseIndex::const_iterator::const_iterator(IndexSlot* data_ptr,
                                            const IndexType& type,
                                            const std::vector<IndexSlot>&
                                            base_vector) :
m_data(data_ptr),
m_type(type),
m_base_begin(base_vector.data()),
m_base_end(base_vector.data()+base_vector.size()) {};

// increment operators
SparseIndex::const_iterator& SparseIndex::const_iterator::operator++() {
    const_iterator& current = *this;
    m_data = find_next_of_type(m_type);
    return current;
};
SparseIndex::const_iterator& SparseIndex::const_iterator::operator++(Integer
                                                                     junk) {
    m_data = find_next_of_type(m_type);
    return *this;
};

// dereference and accessor operators
const Size& SparseIndex::const_iterator::operator*() {
    return (*m_data).first;
};
const Size* SparseIndex::const_iterator::operator->() {
    return &(m_data->first);
};

// comparison operators
bool SparseIndex::const_iterator::operator==(const const_iterator& rhs) {
    return m_data == rhs.m_data; }
bool SparseIndex::const_iterator::operator!=(const const_iterator& rhs) {
    return m_data != rhs.m_data;
};

// static type-finding method
SparseIndex::IndexSlot*
SparseIndex::const_iterator::find_next_of_type(const IndexType& type) const {
    return std::find_if(m_data+1, (IndexSlot*)m_base_end,
                        [&type](IndexSlot& idx_slot){
                            return idx_slot.second == type;
                        });
};
