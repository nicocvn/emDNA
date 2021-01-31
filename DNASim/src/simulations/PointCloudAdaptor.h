// PointCloudAdaptor class
// Nicolas Clauvelin


// point cloud adaptor for nanoflann template implementation

// this class implements the required methods for kd-tree computations

// see PointCloud header for restrictions on CloudType


#ifndef DNASim_PointCloudAdaptor_h
#define DNASim_PointCloudAdaptor_h


#include "simulations/PointCloud.h"


namespace DNASim {


    template <class CloudType> class PointCloudAdaptor {


        // type renaming
        using PointCloud_t = PointCloud<CloudType>;


    public:

        // constructors
        // default constructor is deleted because the const reference requires
        // direct initialization
        PointCloudAdaptor() = delete;
        PointCloudAdaptor(const std::vector<CloudType>& cloud_data);
        PointCloudAdaptor(const PointCloudAdaptor<CloudType>& cloud_adaptor);
        PointCloudAdaptor(PointCloudAdaptor<CloudType>&& cloud_adaptor);

        // copy and move operators
        PointCloudAdaptor<CloudType>&
        operator=(const PointCloudAdaptor<CloudType>& cloud_adaptor);
        PointCloudAdaptor<CloudType>&
        operator=(PointCloudAdaptor<CloudType>&& cloud_adaptor);

        // cloud size accessor
        inline Size kdtree_get_point_count() const {
            return m_cloud.size();
        };

        // distance compuation method
        // p1 is an array with size "size" corresponding to a probe point
        // idx_p2 is the index of the point in the cloud
        // the method has to return the distance between the two points
        // the distance is kept squared for efficiency reasons
        inline Real kdtree_distance(const Real* p1,
                                    const Size idx_p2,
                                    Size size) const {

            // data are assumed to be three dimensional
            DS_ASSERT(size==KdCloudConstants::CloudDim,
                      "wrong dimensionality for cloud data");

            // vectors
            const Vector3 p1_vec({p1[0], p1[2], p1[2]});
            const Vector3 dist_vec = p1_vec-m_cloud.point_coordinates(idx_p2);

            // distance
            return dist_vec.dot(dist_vec);

        };

        // point coordinates accessor
        inline Real kdtree_get_pt(const Size idx, int dim) const {

            // data are assumed to be three dimensional
            DS_ASSERT((0<=dim)&&(dim<=KdCloudConstants::CloudDim),
                      "wrong dimensionality for cloud data");

            return m_cloud.point_coordinate(idx, (Size)dim);

        };

        // bounding-box computation method (defaulted)
        template <class BBOX>
        bool kdtree_get_bbox(BBOX &bb) const {
            return false;
        };


    private:

        PointCloud_t m_cloud;


    };


    // class constructor with data initialization
    template <class CloudType> PointCloudAdaptor<CloudType>::
    PointCloudAdaptor(const std::vector<CloudType>& cloud_data) :
    m_cloud(cloud_data) {};


    // class constructor by copy
    template <class CloudType> PointCloudAdaptor<CloudType>::
    PointCloudAdaptor(const PointCloudAdaptor<CloudType>& cloud_adaptor) :
    m_cloud(cloud_adaptor.m_cloud) {};


    // class constructor by moving
    template <class CloudType> PointCloudAdaptor<CloudType>::
    PointCloudAdaptor(PointCloudAdaptor<CloudType>&& cloud_adaptor) :
    m_cloud(cloud_adaptor.m_cloud) {};


    // copy and move operators
    template <class CloudType>
    PointCloudAdaptor<CloudType>& PointCloudAdaptor<CloudType>::
    operator=(const PointCloudAdaptor<CloudType>& cloud_adaptor) {
        m_cloud = cloud_adaptor.m_cloud;
        return *this;
    };
    template <class CloudType>
    PointCloudAdaptor<CloudType>& PointCloudAdaptor<CloudType>::
    operator=(PointCloudAdaptor<CloudType>&& cloud_adaptor) {
        m_cloud = cloud_adaptor.m_cloud;
        return *this;
    };


}


#endif  // DNASim_PointCloudAdaptor_h
