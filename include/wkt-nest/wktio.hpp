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

  inline size_t read_multiplier(std::istream &in) {
    size_t multiplier;
    while (in)
      switch (in.peek()) {
      case ' ':
        in.get();
        break;
      case 'x':
        in.get();
        in >> multiplier;
        return multiplier;
      default:
        return 1;
      }
    return 1;
  }

  inline
  std::vector<polygon_t> read_polygons(std::istream &in) {
    std::vector<polygon_t> polygons;
    while (in)
      if (auto p = read_polygon(in)) {
        size_t m = read_multiplier(in);
        for (size_t i=0; i<m; i++)
          polygons.push_back(*p);
      }
    return polygons;
  }
}

#endif //WKTNEST_WKTIO_HPP
