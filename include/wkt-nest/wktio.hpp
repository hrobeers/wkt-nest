#pragma once

#ifndef WKTNEST_WKTIO_HPP
#define WKTNEST_WKTIO_HPP

#include <vector>
#include <iostream>
#include <optional>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>

namespace wktnest
{
  typedef boost::geometry::model::d2::point_xy<double> point_type;
  typedef boost::geometry::model::polygon<point_type> polygon;
  typedef boost::geometry::model::box<point_type> box;

  std::optional<box> read_box(std::istream &in);
  std::optional<polygon> read_polygon(std::istream &in);

  inline
  std::vector<polygon> read_polygons(std::istream &in) {
    std::vector<polygon> polygons;
    while (in)
      if (auto p = read_polygon(in))
        polygons.push_back(*p);
    return polygons;
  }
}

#endif //WKTNEST_WKTIO_HPP
