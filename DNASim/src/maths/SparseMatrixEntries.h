// SparseMatrixEntries class
// Nicolas Clauvelin


// wrapper class for vectors of Eigen triplets used to create Eigen sparse
// matrix

// for maximum efficiency the constructor with size initialization should be
// used with a correct estimate of the number of triplets


#ifndef DNASim_SparseMatrixEntries_h
#define DNASim_SparseMatrixEntries_h


#include "maths/Eigen_Includes.h"


namespace DNASim {


    class Matrix3;
    class Matrix4;
    class MatrixN;


    class SparseMatrixEntries {


    public:

        // constructors
        SparseMatrixEntries() = default;
        SparseMatrixEntries(Size n_entries);
        SparseMatrixEntries(const SparseMatrixEntries& entries) = default;
        SparseMatrixEntries(SparseMatrixEntries&& entries) = default;
        ~SparseMatrixEntries() = default;

        // copy and move operators
        SparseMatrixEntries&
        operator=(const SparseMatrixEntries& entries) = default;
        SparseMatrixEntries&
        operator=(SparseMatrixEntries&& entries) = default;

        // add entry methods
        void add_entry(const EigenTriplet& entry);
        void add_entries_from_matrix(const Matrix3& mat3,
                                     const Size& row_start = 0,
                                     const Size& col_start = 0);
        void add_entries_from_matrix(const Matrix4& mat4,
                                     const Size& row_start = 0,
                                     const Size& col_start = 0);
        void add_entries_from_matrix(const MatrixN& mat,
                                     const Size& row_start = 0,
                                     const Size& col_start = 0);

        // clear method
        void clear();

        // iterators
        std::vector<EigenTriplet>::const_iterator begin() const;
        std::vector<EigenTriplet>::const_iterator end() const;


    private:

        // list of triplet
        std::vector<EigenTriplet> m_entries;


    };


}


#endif  // DNASim_SparseMatrixEntries_h
