#pragma once

#ifndef WKTNEST_GEOMETRY_HPP
#define WKTNEST_GEOMETRY_HPP

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

namespace wktnest
{
  typedef double N;
  typedef boost::geometry::model::d2::point_xy<N> point_t;
  typedef boost::geometry::model::polygon<point_t> polygon_t;
  typedef boost::geometry::model::multi_polygon<polygon_t> multi_polygon_t;
  typedef boost::geometry::model::box<point_t> box_t;

  typedef boost::geometry::strategy::transform::matrix_transformer<N,2,2> matrix_t;
  const matrix_t identity_matrix = {1,0,0,
                                    0,1,0,
                                    0,0,1};
}

#endif //WKTNEST_GEOMETRY_HPP
