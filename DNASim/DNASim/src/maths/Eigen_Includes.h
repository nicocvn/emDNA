// Eigen Includes
// Nicolas Clauvelin


// header for Eigen-based implementation


#ifndef DNASim_Eigen_Includes_h
#define DNASim_Eigen_Includes_h


#include <DNASim_Includes.h>


// suppress warnings
#pragma clang diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wshadow"


// includes
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>


// directives
#define EIGEN_NO_DEBUG
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#ifndef EIGEN_FAST_MATH
    #define EIGEN_FAST_MATH
#endif // EIGEN_FAST_MATH

// typedefs
#ifdef DNASIM_SINGLE_PREC
    typedef Eigen::VectorXf EigenVectorX;
    typedef Eigen::Vector4f EigenVector4;
    typedef Eigen::Vector3f EigenVector3;
    typedef Eigen::MatrixXf EigenMatrixX;
    typedef Eigen::Matrix4f EigenMatrix4;
    typedef Eigen::Matrix3f EigenMatrix3;
    typedef Eigen::ArrayXXf EigenArrayXX;
    typedef Eigen::Transform<Real, 3, Eigen::Affine> EigenTransform3D;
    typedef Eigen::Triplet<Real> EigenTriplet;
    typedef Eigen::SparseMatrix<Real> EigenSparseMatrix;
#endif	// DNASIM_SINGLE_PREC
#ifdef DNASIM_DOUBLE_PREC
    typedef Eigen::VectorXd EigenVectorX;
    typedef Eigen::Vector4d EigenVector4;
    typedef Eigen::Vector3d EigenVector3;
    typedef Eigen::MatrixXd EigenMatrixX;
    typedef Eigen::Matrix4d EigenMatrix4;
    typedef Eigen::Matrix3d EigenMatrix3;
    typedef Eigen::ArrayXXd EigenArrayXX;
    typedef Eigen::Transform<Real, 3, Eigen::Affine> EigenTransform3D;
    typedef Eigen::Triplet<Real> EigenTriplet;
    typedef Eigen::SparseMatrix<Real> EigenSparseMatrix;
#endif	// DNASIM_DOUBLE_PREC

// typedefs for eigenvalues and eigenvectors
using Eigenvalues = std::vector<DNASim::RealPair>;
using Eigenvectors = std::vector<std::vector<DNASim::RealPair>>;
using Eigensystem = std::pair<Eigenvalues,Eigenvectors>;
using EigenPair = std::pair<Real,std::vector<Real>>;


#endif  // DNASim_Eigen_Includes_h
