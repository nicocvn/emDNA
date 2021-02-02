// SparseIndex class
// Nicolas Clauvelin


// this class implements a index container with the possibiity to "flag" some
// index values
// const_iterators are provided to iterate over the "clean" and "flagged" index
// values


#ifndef DNASim_SparseIndex_h
#define DNASim_SparseIndex_h


#include "DNASim_Includes.h"


class SparseIndex {


    // index type enum class
    enum class IndexType : bool {
        Clean = true,
        Flagged = false
    };

    // pair type
    using IndexSlot = std::pair<Size,IndexType>;


public:

    // constructors
    // default constructor is deleted to force range initialization
    // range constructor initializes the indices to [start,end) (exclusive)
    SparseIndex() = delete;
    SparseIndex(const Size& index_start, const Size& index_end);
    SparseIndex(const SparseIndex& sparse_index) = default;
    SparseIndex(SparseIndex&& sparse_index) = default;
    ~SparseIndex() = default;

    // copy and move operators
    SparseIndex& operator=(const SparseIndex& sparse_index) = default;
    SparseIndex& operator=(SparseIndex&& sparse_index) = default;

    // index properties
    Size n_index_values() const { return m_indices.size(); };
    Size n_clean_indices() const;
    Size n_flagged_indices() const;
    Size index_first_value() const { return m_indices.front().first; };
    Size index_last_value() const { return m_indices.back().first; };

    // flagging methods
    void flag_index(const Size& index_position) {
        m_indices[index_position].second = IndexType::Flagged;
    };
    void flag_index(const Size& index_start, const Size& index_end) {
        for (Size i=index_start; i<index_end; ++i)
            m_indices[i].second = IndexType::Flagged;
    };

    // forward const_iterator
    class const_iterator {

    public:

        // constructors
        const_iterator() = default;
        const_iterator(IndexSlot* data_ptr,
                       const IndexType& type,
                       const std::vector<IndexSlot>& base_vector);
        const_iterator(const const_iterator& it) = default;
        const_iterator(const_iterator&& it) = default;

        // copy and move operators
        const_iterator& operator=(const const_iterator& it) = default;
        const_iterator& operator=(const_iterator&& it) = default;

        // operators
        const_iterator& operator++();
        const_iterator& operator++(Integer junk);
        const Size& operator*();
        const Size* operator->();

        // comparison operators
        bool operator==(const const_iterator& rhs);
        bool operator!=(const const_iterator& rhs);

    private:

        // static type-finding method
        IndexSlot* find_next_of_type(const IndexType& type) const;
        
        IndexSlot* m_data;
        IndexType m_type;
        const IndexSlot* m_base_begin;
        const IndexSlot* m_base_end;
        
    };
    const_iterator begin_clean() {
        return const_iterator(m_indices.data(),
                              IndexType::Clean,
                              m_indices);
    };
    const_iterator begin_flagged() {
        return const_iterator(m_indices.data(),
                              IndexType::Flagged,
                              m_indices);
    };
    const_iterator end_clean() {
        return const_iterator(m_indices.data()+m_indices.size(),
                              IndexType::Clean,
                              m_indices);
    };
    const_iterator end_flagged() {
        return const_iterator(m_indices.data()+m_indices.size(),
                              IndexType::Flagged,
                              m_indices);
    };


private:

    std::vector<IndexSlot> m_indices;


};


// iterator_traits
namespace std {
template <>
struct iterator_traits<SparseIndex::const_iterator> {
    typedef ptrdiff_t difference_type;
    typedef Size value_type;
    typedef const Size& reference;
    typedef const Size* pointer;
    typedef std::forward_iterator_tag iterator_category;
};
}


#endif  // DNASim_SparseIndex_h
