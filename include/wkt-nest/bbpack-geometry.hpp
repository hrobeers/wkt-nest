#pragma once

#ifndef BBPACK_GEOMETRY_HPP
#define BBPACK_GEOMETRY_HPP

#include "geometry.hpp"

namespace bbpack { namespace geometry {
    using wktnest::matrix_t;
    using wktnest::identity_matrix;

    typedef wktnest::N flt_t;

    typedef int_fast64_t crd_t;
    //typedef flt_t crd_t;

    typedef boost::geometry::model::d2::point_xy<crd_t> point_t;
    typedef boost::geometry::model::polygon<point_t> polygon_t;
    typedef boost::geometry::model::box<point_t> box_t;

    struct dim_t {
      crd_t w;
      crd_t h;
    };

    flt_t area(const box_t& b);
    flt_t area(const polygon_t& b);

    flt_t distance(const point_t& pnt, const polygon_t& poly);

    void centroid(const box_t& b, point_t& c);
    void centroid(const polygon_t& b, point_t& c);

    bool intersects(const box_t& i1, const box_t& i2);
    bool intersects(const box_t& i1, const polygon_t& i2);
    bool intersects(const polygon_t& i1, const polygon_t& i2);

    polygon_t buffer(const polygon_t& p, crd_t buffer_distance);

    inline
    box_t union_(const box_t& b1, const box_t& b2) {
      return {{std::min(b1.min_corner().x(),b2.min_corner().x()),
               std::min(b1.min_corner().y(),b2.min_corner().y())},
              {std::max(b1.max_corner().x(),b2.max_corner().x()),
               std::max(b1.max_corner().y(),b2.max_corner().y())}};
    }
    polygon_t union_(const polygon_t& p1, const polygon_t& p2);
    flt_t union_area(const polygon_t& p1, const polygon_t& p2);
    polygon_t convex_hull(const polygon_t& p1, const polygon_t& p2);

    polygon_t convert(const box_t& b);
}}

#endif //BBPACK_GEOMETRY_HPP
