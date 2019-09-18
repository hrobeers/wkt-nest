#pragma once

#ifndef BBPACK_GEOMETRY_HPP
#define BBPACK_GEOMETRY_HPP

#include "geometry.hpp"

namespace bbpack { namespace geometry {
    using wktnest::matrix_t;
    using wktnest::identity_matrix;

    typedef wktnest::N flt_t;

    typedef int_fast32_t crd_t;
    //typedef flt_t crd_t;

    typedef boost::geometry::model::d2::point_xy<crd_t> point_t;
    typedef boost::geometry::model::polygon<point_t> polygon_t;
    typedef boost::geometry::model::box<point_t> box_t;

    struct dim_t {
      crd_t w;
      crd_t h;
    };

    flt_t area(const box_t& b);

    void centroid(const box_t& b, point_t& c);
    void centroid(const polygon_t& b, point_t& c);

    bool intersects(const box_t& i1, const box_t& i2);
    bool intersects(const box_t& i1, const polygon_t& i2);
    bool intersects(const polygon_t& i1, const polygon_t& i2);
}}

#endif //BBPACK_GEOMETRY_HPP
