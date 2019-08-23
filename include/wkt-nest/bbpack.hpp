#pragma once

#ifndef WKTNEST_BBPACK_HPP
#define WKTNEST_BBPACK_HPP

#include <vector>
#include <map>

#include "wkt-nest/geometry.hpp"

namespace wktnest {
  namespace bbpack {
    struct node_t {
      box_t box;
      bool used = false;
      node_t* up = nullptr;
      node_t* right = nullptr;
    };
    struct state {
      box_t bin;
      // list of nodes per bin
      //std::vector<std::vector<node_t>> nodes;
      // single bin of nodes TODO switch to multiple bins
      std::vector<node_t> nodes;
      // fit per bbox
      std::map<const box_t*,node_t*> fits;
    };

    state init(box_t& bin);
    std::vector<box_t> fit(state& s, std::vector<polygon_t> polygons);
    node_t* find_node(state& s, node_t* root, box_t& bbox);
    node_t* split_node(state& s, node_t* n, box_t& bbox);
  }
}

#endif //WKTNEST_BBPACK_HPP
