#pragma once

#ifndef WKTNEST_GEOMETRY_HPP
#define WKTNEST_GEOMETRY_HPP

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace wktnest
{
  typedef double flt_t;
  typedef int_fast32_t crd_t;

  typedef boost::geometry::model::d2::point_xy<crd_t> point_t;
  typedef boost::geometry::model::polygon<point_t> polygon_t;
  typedef boost::geometry::model::box<point_t> box_t;

  namespace ublas = boost::numeric::ublas;
  typedef std::vector<flt_t> array_t;
  typedef ublas::matrix<flt_t, ublas::row_major, array_t> matrix_t;
  typedef ublas::vector<flt_t, array_t> vector_t;
  typedef ublas::identity_matrix<flt_t> identity_matrix;
}

#endif //WKTNEST_GEOMETRY_HPP
