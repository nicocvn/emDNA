// BpTwistDensity class
// Nicolas Clauvelin


#include <Segment.h>
#include <Triad.h>
#include <CurveTopology.h>
#include <BpStepTwistDensity.h>


namespace DNASim {


    // static methods for twist density computation
    // the first and last values are zero as the computation requires
    // "internal" data (i.e., the first and last vertex are incomplete)
    std::vector<Vector3>
    BpTwistDensity::compute_twist_density(const std::vector<BasePair>&
                                          base_pairs) {

        // number of base pairs
        Size n_bp = base_pairs.size();

        // collect the base-pair origins to form the centerline
        std::vector<Vertex> centerline;
        centerline.reserve(n_bp);
        for (Size i=0; i<n_bp; ++i)
            centerline.push_back(base_pairs[i].origin());
        const std::vector<Segment>
        centerline_segments(CurveTopology::curve_segments(centerline));

        // vertex frames
        // this frame is define by the vertex tangent and the binormal
        std::vector<Triad> vertex_frames(n_bp);
        vertex_frames[0] = compute_vertex_frame(centerline_segments[n_bp-1],
                                                centerline_segments[0]);
        for (Size i=0; i<n_bp-2; ++i)
            vertex_frames[i+1] = compute_vertex_frame(centerline_segments[i],
                                                      centerline_segments[i+1]);
        vertex_frames[n_bp-1] =
        compute_vertex_frame(centerline_segments[n_bp-2],
                             centerline_segments[n_bp-1]);

        // twist densities
        std::vector<Vector3> twist_densities;
        twist_densities.reserve(n_bp);
        for (Size i=0; i<n_bp-1; ++i) {
            twist_densities.
            push_back(compute_bp_step_twist_density(base_pairs[i],
                                                    base_pairs[i+1],
                                                    vertex_frames[i],
                                                    vertex_frames[i+1]));
        };
        twist_densities.
        push_back(compute_bp_step_twist_density(base_pairs[n_bp-1],
                                                base_pairs[0],
                                                vertex_frames[n_bp-1],
                                                vertex_frames[0]));

        return twist_densities;

    };


    // bp step local twist density
    Vector3
    BpTwistDensity::compute_bp_step_twist_density(const BasePair& bp1,
                                                  const BasePair& bp2,
                                                  const Triad& vertex_frame1,
                                                  const Triad& vertex_frame2) {

        // projection of the base-pair frames axes onto the vertex frames axes
        Vector3 proj_bp1_I(bp1.axis(I).dot(vertex_frame1.axis(I))*
                           vertex_frame1.axis(I)
                           +
                           bp1.axis(I).dot(vertex_frame1.axis(J))*
                           vertex_frame1.axis(J));
        proj_bp1_I.normalize();
        Vector3 proj_bp1_J(bp1.axis(J).dot(vertex_frame1.axis(I))*
                           vertex_frame1.axis(I)
                           +
                           bp1.axis(J).dot(vertex_frame1.axis(J))*
                           vertex_frame1.axis(J));
        proj_bp1_J.normalize();
        Vector3 proj_bp2_I(bp2.axis(I).dot(vertex_frame2.axis(I))*
                           vertex_frame2.axis(I)
                           +
                           bp2.axis(I).dot(vertex_frame2.axis(J))*
                           vertex_frame2.axis(J));
        proj_bp2_I.normalize();
        Vector3 proj_bp2_J(bp2.axis(J).dot(vertex_frame2.axis(I))*
                           vertex_frame2.axis(I)
                           +
                           bp2.axis(J).dot(vertex_frame2.axis(J))*
                           vertex_frame2.axis(J));
        proj_bp2_J.normalize();

        // projection base-pair frames
        Triad proj_bp1;
        proj_bp1.set_origin(vertex_frame1.origin());
        proj_bp1.set_axis(I, proj_bp1_I);
        proj_bp1.set_axis(J, proj_bp1_J);
        proj_bp1.set_axis(K, vertex_frame1.axis(K));
        Triad proj_bp2;
        proj_bp2.set_origin(vertex_frame2.origin());
        proj_bp2.set_axis(I, proj_bp2_I);
        proj_bp2.set_axis(J, proj_bp2_J);
        proj_bp2.set_axis(K, vertex_frame2.axis(K));

        // binormal angle
        const Vector3 step_tangent(Segment(bp1.origin(),
                                           bp2.origin()).tangent());
        const Real cosval(vertex_frame1.axis(J).
                          dot(vertex_frame2.axis(J)));
        const Real sinval(vertex_frame1.axis(J).
                          cross(vertex_frame2.axis(J)).dot(step_tangent));
        const Real binormal_angle(std::atan2(sinval, cosval));

        // axes
        const Vector3 binormal1 = vertex_frame1.axis(J);
        const Vector3 binormal2 = vertex_frame2.axis(J);
        const Vector3 tangent1 = vertex_frame1.axis(K);
        const Vector3 tangent2 = vertex_frame2.axis(K);

        // short axis twist density
        const Real short_sinval1 =
        binormal1.cross(proj_bp1.axis(I)).dot(tangent1);
        const Real short_cosval1 = binormal1.dot(proj_bp1.axis(I));
        const Real short_ang1 = std::atan2(short_sinval1, short_cosval1);
        const Real short_sinval2 =
        binormal2.cross(proj_bp2.axis(I)).dot(tangent2);
        const Real short_cosval2 = binormal2.dot(proj_bp2.axis(I));
        const Real short_ang2 = std::atan2(short_sinval2, short_cosval2);
        Real short_tw = short_ang2 - short_ang1 + binormal_angle;

        // long axis twist density
        const Real long_sinval1 =
        binormal1.cross(proj_bp1.axis(J)).dot(tangent1);
        const Real long_cosval1 = binormal1.dot(proj_bp1.axis(J));
        const Real long_ang1 = std::atan2(long_sinval1, long_cosval1);
        const Real long_sinval2 =
        binormal2.cross(proj_bp2.axis(J)).dot(tangent2);
        const Real long_cosval2 = binormal2.dot(proj_bp2.axis(J));
        const Real long_ang2 = std::atan2(long_sinval2, long_cosval2);
        Real long_tw = long_ang2 - long_ang1 + binormal_angle;

        // average axis twist density
        Vector3 avgaxis1(proj_bp1.axis(I)+proj_bp1.axis(J));
        avgaxis1.normalize();
        Vector3 avgaxis2(proj_bp2.axis(I)+proj_bp2.axis(J));
        avgaxis2.normalize();
        const Real avg_sinval1 = binormal1.cross(avgaxis1).dot(tangent1);
        const Real avg_cosval1 = binormal1.dot(avgaxis1);
        const Real avg_ang1 = std::atan2(avg_sinval1, avg_cosval1);
        const Real avg_sinval2 = binormal2.cross(avgaxis2).dot(tangent2);
        const Real avg_cosval2 = binormal2.dot(avgaxis2);
        const Real avg_ang2 = std::atan2(avg_sinval2, avg_cosval2);
        Real avg_tw = avg_ang2 - avg_ang1 + binormal_angle;

        // twist correction
        if (short_tw > F_PI)
            short_tw -= Real(2)*F_PI;
        else if (short_tw < -F_PI)
            short_tw += Real(2)*F_PI;
        if (long_tw > F_PI)
            long_tw -= Real(2)*F_PI;
        else if (long_tw < -F_PI)
            long_tw += Real(2)*F_PI;
        if (avg_tw > F_PI)
            avg_tw -= Real(2)*F_PI;
        else if (avg_tw < -F_PI)
            avg_tw += Real(2)*F_PI;

        return Vector3(short_tw, long_tw, avg_tw);

    };


    // vertex frame computation method
    Triad BpTwistDensity::compute_vertex_frame(const Segment& ingoing,
                                               const Segment& outgoing) {

        // vertex binormal
        Vector3 binormal(ingoing.tangent().cross(outgoing.tangent()));
        binormal.normalize();

        // vertex tangent
        Vector3 tangent(ingoing.tangent()+outgoing.tangent());
        tangent.normalize();

        // vertex normal
        Vector3 normal(binormal.cross(tangent));
        normal.normalize();

        // vertex frame
        Triad vertex_frame;
        vertex_frame.set_origin(ingoing.second());
        vertex_frame.set_axis(I, normal);
        vertex_frame.set_axis(J, binormal);
        vertex_frame.set_axis(K, tangent);

        return vertex_frame;

    };


}

