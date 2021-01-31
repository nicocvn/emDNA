// PointCloudKdTree class
// Nicolas Clauvelin


// kd-tree template class implementation for point cloud data

// the max leaf setting is set to its default value

// see PointCloud header for restrictions on CloudType

// caution:
// all distances returned from kd methods are squared!!!


#ifndef DNASim_PointCloudKdTree_h
#define DNASim_PointCloudKdTree_h


#include "simulations/PointCloudAdaptor.h"
#include <nanoflann.hpp>


namespace DNASim {


    // types renaming
    // the Size data refers to the indices of points in the point cloud
    using NeighborResults = std::tuple<std::vector<Size>,std::vector<Real>>;
    using RadiusResults = std::vector<std::pair<Size,Real>>;


    template <class CloudType> class PointCloudKdTree {


        // type renaming
        // use L2_Adaptor cause L2_Simple causes a bug in distance computation
        using PointCloudAdaptor_t = PointCloudAdaptor<CloudType>;
        using CloudKdTree_t =
        nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Adaptor<Real, PointCloudAdaptor_t>,
            PointCloudAdaptor_t,
            KdCloudConstants::CloudDim
        >;


    public:

        // constructors
        // default constructor is deleted because the const reference requires
        // direct initialization
        PointCloudKdTree() = delete;
        PointCloudKdTree(const std::vector<CloudType>& data);

        // kd-tree method - nearest neighbors
        NeighborResults
        find_nearest_neighbors_from_point(const Vector3& point,
                                          const Size& n_neighbors,
                                          const bool& sorted_search) const;

        // kd-tree method - radius search
        // if n_points_estimate is set to zero no reserve is peformed
        RadiusResults
        radius_search_from_point(const Vector3& point,
                                 const Real& radius,
                                 const Size& n_points_estimate,
                                 const bool& sorted_search) const;


    private:

        const PointCloudAdaptor_t m_cloud_adaptor;
        CloudKdTree_t m_kd_tree;


    };


    // class constructor with data initialization
    template <class CloudType> PointCloudKdTree<CloudType>::
    PointCloudKdTree(const std::vector<CloudType>& cloud_data) :
    m_cloud_adaptor(cloud_data),
    m_kd_tree(KdCloudConstants::CloudDim, m_cloud_adaptor) {
        m_kd_tree.buildIndex();
    };


    // nearest neighbors method
    template <class CloudType>
    NeighborResults PointCloudKdTree<CloudType>::
    find_nearest_neighbors_from_point(const Vector3& point,
                                      const Size& n_neighbors,
                                      const bool& sorted_search) const {

        // results set and result containers
        nanoflann::KNNResultSet<Real> results_set(n_neighbors);
        std::vector<Size> neighbor_indices(n_neighbors, Size(0));
        std::vector<Real> neighbor_sq_distances(n_neighbors, FLOAT_INIT);
        results_set.init(neighbor_indices.data(),
                         neighbor_sq_distances.data());

        // query point
        const Real query_point[KdCloudConstants::CloudDim] =
        {point[X], point[Y], point[Z]};

        // search params
        nanoflann::SearchParams search_params;
        search_params.checks = 0; // not used so whatever
        search_params.eps = 0; // exact search
        search_params.sorted = sorted_search;

        // search
        m_kd_tree.findNeighbors(results_set,
                                &query_point[0],
                                search_params);

        return NeighborResults(neighbor_indices, neighbor_sq_distances);

    };


    // radius search method
    template <class CloudType>
    RadiusResults PointCloudKdTree<CloudType>::
    radius_search_from_point(const Vector3& point,
                             const Real& radius,
                             const Size& n_points_estimate,
                             const bool& sorted_search) const {

        // radius squared
        const Real radius_squared = radius*radius;

        // query point
        const Real query_point[KdCloudConstants::CloudDim] =
        {point[X], point[Y], point[Z]};

        // container
        std::vector<std::pair<Size,Real>> results;
        if (n_points_estimate != 0)
            results.reserve(n_points_estimate);

        // search params
        nanoflann::SearchParams search_params;
        search_params.checks = 0; // not used so whatever
        search_params.eps = 0; // exact search
        search_params.sorted = sorted_search;

        // search
        m_kd_tree.radiusSearch(&query_point[0], radius_squared,
                               results,
                               search_params);

        return results;

    };


};


#endif  // DNASim_PointCloudKdTree_h
