#pragma once

#ifndef WKTNEST_GEOMETRY_HPP
#define WKTNEST_GEOMETRY_HPP

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace wktnest
{
  typedef boost::geometry::model::d2::point_xy<double> point_t;
  typedef boost::geometry::model::polygon<point_t> polygon_t;
  typedef boost::geometry::model::box<point_t> box_t;
}

#endif //WKTNEST_GEOMETRY_HPP
