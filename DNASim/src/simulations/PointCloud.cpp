// PointCloud class
// Nicolas Clauvelin


// the following methods correspond to specializations of the coordinates
// accessors


#include <Shape.h>
#include <PointCloud.h>


namespace DNASim {


    // point coordinate accessor - Vector3 specialization
    template <>
    Real PointCloud<Vector3>::point_coordinate(const Size& idx,
                                               const Size& coord) const {
        return m_cloud_points[idx][static_cast<CartesianCoordinates>(coord)];
    };


    // point coordinate accessor - Triad specialization
    template <>
    Real PointCloud<Triad>::point_coordinate(const Size& idx,
                                             const Size& coord) const {
        const Vector3& origin = m_cloud_points[idx].origin();
        return origin[static_cast<CartesianCoordinates>(coord)];
    };


    // point coordinate accessor - Shape_Ptr specialization
    template <>
    Real PointCloud<Shape_Ptr>::point_coordinate(const Size& idx,
                                                 const Size& coord) const {
        const Vector3& origin = m_cloud_points[idx]->position();
        return origin[static_cast<CartesianCoordinates>(coord)];
    };


    // three-dimensional coordinates accessor - Vector3 specialization
    template <>
    Vector3 PointCloud<Vector3>::point_coordinates(const Size& idx) const
    {
        return m_cloud_points[idx];
    };


    // three-dimensional coordinates accessor - Triad specialization
    template <>
    Vector3 PointCloud<Triad>::point_coordinates(const Size& idx) const {
        return m_cloud_points[idx].origin();
    };


    // three-dimensional coordinates accessor - Shape_Ptr specialization
    template <>
    Vector3 PointCloud<Shape_Ptr>::point_coordinates(const Size& idx) const {
        return m_cloud_points[idx]->position();
    };


}
