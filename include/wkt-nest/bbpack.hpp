#pragma once

#ifndef WKTNEST_BBPACK_HPP
#define WKTNEST_BBPACK_HPP

#include <vector>
#include <list>
#include <map>

#include "wkt-nest/geometry.hpp"

namespace wktnest {
  // TODO move to own header when reused elsewhere
  enum SORTING { NONE, HEIGHT, SHUFFLE };
  struct nesting_opts {
    bool bbox = false;
    SORTING sorting = NONE;
  };

  namespace bbpack {
    struct placement {
      size_t bin= 0;
      polygon_t polygon;
      box_t bbox;
      matrix_t transform;
    };
    typedef std::vector<placement> fit_result;

    fit_result fit(const box_t& bin, const std::vector<polygon_t>& polygons, SORTING sorting, bool compact);

    template<typename T_opts>
    fit_result fit(const box_t& bin, const std::vector<polygon_t>& polygons, const T_opts& opts) {
      return fit(bin, polygons, opts.sorting, !opts.bbox);
    }
  }
}

#endif //WKTNEST_BBPACK_HPP
