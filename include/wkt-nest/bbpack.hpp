#pragma once

#ifndef WKTNEST_BBPACK_HPP
#define WKTNEST_BBPACK_HPP

#include <vector>
#include <list>
#include <map>

#include "wkt-nest/geometry.hpp"

namespace wktnest {
  // TODO move to own header when reused elsewhere
  enum SORTING { NONE, HEIGHT };
  struct nesting_opts {
    SORTING sorting = NONE;
  };

  namespace bbpack {
    class item_t {
      bool _placed = false;
      const polygon_t* _source;
      matrix_t _init_transform;

      polygon_t _polygon;
      box_t _bbox;
      matrix_t _transform;

    public:
      item_t(const polygon_t* p);

      const polygon_t* source() const { return _source; }
      const polygon_t* polygon() const { return &_polygon; }
      const box_t* bbox() const { return &_bbox; }
      const matrix_t* transform() const { return &_transform; }
      void placed(bool p) { _placed = p; }
      bool placed() const { return _placed; }

      void relative_transform(const matrix_t& t);
      void absolute_transform(const matrix_t& t);
    };
    struct node_t {
      box_t box;
      bool used = false;
      node_t* up = nullptr;
      node_t* right = nullptr;
    };
    struct state_t {
      box_t bin;
      SORTING sorting;
      std::vector<item_t> items;
      std::map<const polygon_t*, const item_t*> fits;
      std::map<const item_t*, size_t> item_to_bin_idx;

      // TODO below still needed?

      // list of nodes per bin
      //std::vector<std::list<node_t>> nodes;
      // single bin of nodes TODO switch to multiple bins?
      // Nodes should be stored in list since vector will reallocate on resizing
      std::list<node_t> nodes;
    };

    template<typename T_opts>
    state_t init(const box_t& bin, const T_opts& opts) {
      return { bin, opts.sorting };
    }

    std::vector<matrix_t> fit(state_t& s, const std::vector<polygon_t>& polygons);
    node_t* find_node(state_t& s, node_t* root, item_t* item);
    node_t* split_node(state_t& s, node_t* n, item_t* item);
  }
}

#endif //WKTNEST_BBPACK_HPP
