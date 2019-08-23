#pragma once

#ifndef WKTNEST_WKTIO_HPP
#define WKTNEST_WKTIO_HPP

#include <vector>
#include <iostream>
#include <optional>

#include "wkt-nest/geometry.hpp"

namespace wktnest
{
  std::optional<box_t> read_box(std::istream &in);
  std::optional<polygon_t> read_polygon(std::istream &in);

  inline
  std::vector<polygon_t> read_polygons(std::istream &in) {
    std::vector<polygon_t> polygons;
    while (in)
      if (auto p = read_polygon(in))
        polygons.push_back(*p);
    return polygons;
  }
}

#endif //WKTNEST_WKTIO_HPP
