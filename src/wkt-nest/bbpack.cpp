#include "wkt-nest/bbpack.hpp"

#include <algorithm>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/geometry/algorithms/for_each.hpp>

namespace bg = boost::geometry;

using namespace wktnest;
using namespace bbpack;

namespace {
  struct dim_t {
    double w;
    double h;
  };
  dim_t dims(const box_t* b) {
    double min_x = bg::get<bg::min_corner, 0>(*b);
    double min_y = bg::get<bg::min_corner, 1>(*b);
    double max_x = bg::get<bg::max_corner, 0>(*b);
    double max_y = bg::get<bg::max_corner, 1>(*b);
    return { max_x-min_x, max_y-min_y };
  }
}


bbpack::item_t::item_t(const polygon_t& p) {
  _source = &p;
  // TODO apply transform identity instead of code below
  // TODO transform to origin
  _polygon = p;
  bg::envelope(_polygon,_bbox);
}


state_t bbpack::init(const box_t& bin) {
  return { bin };
}

std::vector<box_t> bbpack::fit(state_t& s, const std::vector<polygon_t>& polygons) {

  for (auto p : polygons)
    s.items.push_back(item_t(p));

  // construct root
  s.nodes.push_back({s.bin});
  node_t* root = &s.nodes.back();

  for (item_t& item : s.items)
    if (auto node = find_node(s, root, &item)) // TODO else pack in other bin
      s.fits[&item] = split_node(s, node, &item);

  // Create fitted box vector
  std::vector<box_t> fboxes;
  fboxes.resize(s.items.size());
  std::transform(s.items.cbegin(), s.items.cend(), fboxes.begin(),
                 [&s](const item_t& item) -> box_t {
                   // TODO filter first on fit
                   node_t* n = s.fits[&item];
                   //if (!n) return box_t();
                   point_t min = item.bbox()->min_corner();
                   point_t max = item.bbox()->max_corner();
                   bg::add_point(min, n->box.min_corner());
                   bg::add_point(max, n->box.min_corner());
                   return box_t(min, max);
                 });

  return fboxes;
}

node_t* bbpack::find_node(state_t& s, node_t* root, item_t* item) {

  if (root->used) {
    node_t* up = root->up;
    node_t* right = root->right;
    if (auto node = find_node(s, right, item))
      return node;
    return find_node(s, up, item);
  }

  auto rtdims = dims(&root->box);
  auto bbdims = dims(item->bbox());

  if ((bbdims.w <= rtdims.w) && (bbdims.h <= rtdims.h))
    return root;

  return nullptr;
}

node_t* bbpack::split_node(state_t& s, node_t* n, item_t* item) {

  s.fits[item] = n;
  n->used = true;

  double n_min_x = bg::get<bg::min_corner, 0>(n->box);
  double n_min_y = bg::get<bg::min_corner, 1>(n->box);
  double n_max_x = bg::get<bg::max_corner, 0>(n->box);
  double n_max_y = bg::get<bg::max_corner, 1>(n->box);

  auto bbdims = dims(item->bbox());

  // node up
  box_t up = {{n_min_x, n_min_y + bbdims.h}, {n_max_x, n_max_y}};
  s.nodes.push_back({ up });
  n->up = &s.nodes.back();

  // node right
  box_t right = {{n_min_x + bbdims.w, n_min_y}, {n_max_x, n_max_y}};
  s.nodes.push_back({ right });
  n->right = &s.nodes.back();

  return n;
}
