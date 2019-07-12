// PointCloud class
// Nicolas Clauvelin


// point cloud template implementation for three-dimensional data

// the cloud can used arbitray data with the following conditions:
// - the data "positions" have to be three dimensional,
// - the methods point_coordinate/point_coordinates need to be specialized.

// the data are stored as a const reference (and, hence, are not copied)


#ifndef DNASim_PointCloud_h
#define DNASim_PointCloud_h


#include <Vector3.h>
#include <VectorN.h>
#include <Matrix3.h>
#include <MatrixN.h>
#include <Triad.h>


namespace DNASim {


    // internal constants
    namespace KdCloudConstants {
        const Size CloudDim = 3;
    };


    // point cloud interface
    template <class CloudType> class PointCloud {


        // types renaming
        using CloudType_vec = std::vector<CloudType>;


    public:

        // constructors
        // default constructor is deleted because the const reference requires
        // direct initialization
        PointCloud() = delete;
        PointCloud(const CloudType_vec& cloud_data);
        PointCloud(const PointCloud<CloudType>& point_cloud);
        PointCloud(PointCloud<CloudType>&& point_cloud);

        // copy and move operators
        PointCloud<CloudType>&
        operator=(const PointCloud<CloudType>& point_cloud);
        PointCloud<CloudType>&
        operator=(PointCloud<CloudType>&& point_cloud);

        // cloud size accessor
        inline Size size() const {
            return m_cloud_points.size();
        };

        // point coordinates accessor
        // return the coord-th coordinate of the idx-th point
        Real point_coordinate(const Size& idx, const Size& coord) const;

        // three-dimensional coordinates accessor
        Vector3 point_coordinates(const Size& idx) const;

        // data accessors
        const CloudType& point(const Size& idx) const;

        // PC analysis methods
        Vector3 cloud_center() const;
        MatrixN cloud_points_matrix() const;
        Matrix3 cloud_points_covariance_matrix() const;
        Triad cloud_pca_frame() const;


    private:

        const CloudType_vec& m_cloud_points;


    };


    // class constructor with data initialization
    template <class CloudType>
    PointCloud<CloudType>::PointCloud(const CloudType_vec& cloud_data) :
    m_cloud_points(cloud_data) {};


    // class constructor by copy
    template <class CloudType>
    PointCloud<CloudType>::PointCloud(const PointCloud<CloudType>&
                                      point_cloud) :
    m_cloud_points(point_cloud.m_cloud_points) {};


    // class constructor by moving
    template <class CloudType>
    PointCloud<CloudType>::PointCloud(PointCloud<CloudType>&& point_cloud) :
    m_cloud_points(point_cloud.m_cloud_points) {};


    // copy and move operators
    template <class CloudType>
    PointCloud<CloudType>& PointCloud<CloudType>::
    operator=(const PointCloud<CloudType>& point_cloud) {
        m_cloud_points = point_cloud.m_cloud_points;
        return *this;
    };
    template <class CloudType>
    PointCloud<CloudType>& PointCloud<CloudType>::
    operator=(PointCloud<CloudType>&& point_cloud) {
        m_cloud_points = point_cloud.m_cloud_points;
        return *this;
    };


    // data accessor
    template <class CloudType>
    const CloudType& PointCloud<CloudType>::point(const Size& idx) const {
        return m_cloud_points[idx];
    };


    // PC analysis methods - cloud center
    template <class CloudType>
    Vector3 PointCloud<CloudType>::cloud_center() const {

        // sizing
        const Size n_pts = m_cloud_points.size();

        // coordinates accumulation
        Vector3 accumulated_coords;
        for (Size i=0; i<n_pts; ++i)
            accumulated_coords += point_coordinates(i);

        // mean
        return accumulated_coords*(1.0/n_pts);

    };


    // PC analysis methods - cloud points matrix
    template <class CloudType>
    MatrixN PointCloud<CloudType>::cloud_points_matrix() const {

        // sizing
        const Size n_pts = m_cloud_points.size();

        // matrix storage
        MatrixN mat(KdCloudConstants::CloudDim,n_pts);
        for (Size i=0; i<n_pts; ++i) {
            const Vector3& pt = point_coordinates(i);
            mat.set_col(i,
                        VectorN({pt[X], pt[Y], pt[Z]})
                        );
        };

        return mat;

    };


    // PC analysis methods - cloud points covariance matrix
    template <class CloudType>
    Matrix3 PointCloud<CloudType>::cloud_points_covariance_matrix() const {

        // sizing
        const Size n_pts = m_cloud_points.size();

        // geometric center
        const Vector3 center = cloud_center();
        MatrixN center_mat(KdCloudConstants::CloudDim,n_pts);
        for (Size i=0; i<n_pts; ++i)
            center_mat.set_col(i,
                               VectorN({center[X],center[Y],center[Z]}));

        // cloud matrices
        MatrixN mat = cloud_points_matrix()-center_mat;
        MatrixN mat_T = cloud_points_matrix()-center_mat;
        mat_T.transpose();
        mat *= mat_T;

        // size checking
        DS_ASSERT(mat.n_rows()==KdCloudConstants::CloudDim
                  &&
                  mat.n_cols()==KdCloudConstants::CloudDim,
                  "wrong sizing for cloud covariance matrix");

        // convert to Matrix3
        Matrix3 cov_mat;
        for (Size i=0; i<KdCloudConstants::CloudDim; ++i)
            for (Size j=0; j<KdCloudConstants::CloudDim; ++j)
                cov_mat(i,j) = mat(i,j);

        return static_cast<Real>(1.0/(n_pts-1.0))*cov_mat;

    };


    // PC analysis methods - cloud PCA frame
    template <class CloudType>
    Triad PointCloud<CloudType>::cloud_pca_frame() const {

        // covariance matrix
        const Matrix3 cov_mat = cloud_points_covariance_matrix();

        // eigensystem
        Eigensystem eigen_system = cov_mat.eigensystem();

        // the covariance matrix is assumed to be symmetric, so the eigenvalues
        // and eigenvectors have null imaginary parts
        //
        // we collect the eigenvalues and eigenvectors checking for this
        // assumption
        Vector3 eigen_values;
        std::vector<Vector3> eigen_vectors;
        for (Size i=0; i<KdCloudConstants::CloudDim; ++i) {

            // value
            const RealPair val = (eigen_system.first)[i];

            // vector
            const std::vector<RealPair> vec = (eigen_system.second)[i];

            // checking
            DS_ASSERT(std::fabs(val.second)<ZERO_EPS,
                      "non-zero imaginary part in eigenvalues of covariance "
                      "matrix");
            for (Size k=0; k<KdCloudConstants::CloudDim; ++k)
                DS_ASSERT(std::fabs(vec[k].second)<ZERO_EPS,
                          "non-zero imaginary part in eigenvectors of "
                          "covariance matrix");

            // storage
            eigen_values[static_cast<CartesianCoordinates>(i)] = val.first;
            eigen_vectors.push_back(Vector3({
                                            vec[0].first,
                                            vec[1].first,
                                            vec[2].first
                                            })
                                    );

        };

        // reorientation step
        Vector3 u1 = eigen_vectors[0];
        Vector3 u2 = eigen_vectors[1];
        if (point_coordinates(0).dot(u1) < 0)
            u1 *= Real(-1.0);
        if (point_coordinates(0).dot(u2) < 0)
            u2 *= Real(-1.0);
        Vector3 u3 = u1.cross(u2);

        // PCA frame
        Triad pca_frame;
        pca_frame.set_origin(cloud_center());
        pca_frame.set_axis(I, u1);
        pca_frame.set_axis(J, u2);
        pca_frame.set_axis(K, u3);
        pca_frame.orthogonalize();

        return pca_frame;

    };


}


#endif  // DNASim_PointCloud_h
