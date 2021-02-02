// test_PointCloudKdTree class
// Nicolas Clauvelin


#include "test_PointCloudKdTree.h"
#include <DNASim.h>
using namespace DNASim;


namespace {


    // GlobalTest
    TEST_F(PointCloudKdTreeTest, GlobalTest) {

        // type renaming
        using RandomReal = PRN_Real<RandomSeed::No>;

        // random cloud generation
        std::vector<Vector3> cloud;
        cloud.reserve(1000);
        for (Size i=0; i<1000; ++i) {
            Vector3 point(RandomReal::generate_with_bounds(Real(-10.0),
                                                           Real(10.0)),
                          RandomReal::generate_with_bounds(Real(-10.0),
                                                           Real(10.0)),
                          RandomReal::generate_with_bounds(Real(-10.0),
                                                           Real(10.0)));
            cloud.push_back(point);
        };

        DS_ASSERT(cloud.size()==1000, "cloud size error");

        // kd tree
        PointCloudKdTree<Vector3> cloud_kd_tree(cloud);

        // query point
        Vector3 point(RandomReal::generate_with_bounds(Real(-10.0),
                                                       Real(10.0)),
                      RandomReal::generate_with_bounds(Real(-10.0),
                                                       Real(10.0)),
                      RandomReal::generate_with_bounds(Real(-10.0),
                                                       Real(10.0)));

        // nearest neighbors search
        NeighborResults nn_res =
        cloud_kd_tree.find_nearest_neighbors_from_point(point, 10, true);

        // radius search
        RadiusResults rad_res =
        cloud_kd_tree.radius_search_from_point(point, Real(10.0), 0, true);

        // neighbors indices and distances
        std::vector<Size> n_indices = std::get<0>(nn_res);
        std::vector<Real> n_distances = std::get<1>(nn_res);

        // comparison
        // both search are sorted
        for (Size i=0; i<10; ++i) {

            EXPECT_TRUE(rad_res[i].second <= Real(10.0)*Real(10.0));
            EXPECT_EQ(n_indices[i], rad_res[i].first);
            EXPECT_EQ(n_distances[i], rad_res[i].second);

        };

        // distances checking
        for (Size i=0; i<10; ++i) {

            // real distance
            const Real sq_dist =
            (point-cloud[rad_res[i].first]).dot(point-cloud[rad_res[i].first]);

            EXPECT_NEAR(sq_dist, rad_res[i].second, ZERO_EPS);

        };

    };


    // PCAFrameTest
    TEST_F(PointCloudKdTreeTest, PCAFrameTest) {

        // type renaming
        using RandomReal = PRN_Real<RandomSeed::No>;

        // random cloud generation
        std::vector<Vector3> cloud;
        cloud.reserve(10);
        for (Size i=0; i<10; ++i) {
            Vector3 point(RandomReal::generate_with_bounds(Real(-1.0),
                                                           Real(1.0)),
                          RandomReal::generate_with_bounds(Real(-1.0),
                                                           Real(1.0)),
                          RandomReal::generate_with_bounds(Real(-1.0),
                                                           Real(1.0)));
            cloud.push_back(point);
        };

        DS_ASSERT(cloud.size()==10, "cloud size error");

        // pca frame
        PointCloud<Vector3> point_cloud(cloud);
        const Triad pca_frame = point_cloud.cloud_pca_frame();

        // transformation
        std::vector<Vector3> aff_t_cloud;
        aff_t_cloud.reserve(cloud.size());
        for (const Vector3& pt : cloud)
            aff_t_cloud.push_back(pca_frame.express_point_in(pt));

        // compute new pca frame for checking
        PointCloud<Vector3> trans_cloud(aff_t_cloud);
        const Triad aff_t_pca = trans_cloud.cloud_pca_frame();

        // checking
        EXPECT_NEAR(FLOAT_INIT,
                    aff_t_pca.origin().norm(),
                    ZERO_EPS);
        EXPECT_NEAR(Real(1.0),
                    aff_t_pca.axis(I)[X],
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    aff_t_pca.axis(I)[Y],
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    aff_t_pca.axis(I)[Z],
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    aff_t_pca.axis(J)[X],
                    ZERO_EPS);
        EXPECT_NEAR(Real(1.0),
                    aff_t_pca.axis(J)[Y],
                    ZERO_EPS);
        EXPECT_NEAR(FLOAT_INIT,
                    aff_t_pca.axis(J)[Z],
                    ZERO_EPS);

    };


}

