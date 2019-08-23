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
  dim_t dims(const box_t& b) {
    double min_x = bg::get<bg::min_corner, 0>(b);
    double min_y = bg::get<bg::min_corner, 1>(b);
    double max_x = bg::get<bg::max_corner, 0>(b);
    double max_y = bg::get<bg::max_corner, 1>(b);
    return { max_x-min_x, max_y-min_y };
  }
}

state bbpack::init(box_t& bin) {
  return { bin };
}

std::vector<box_t> bbpack::fit(state& s, std::vector<polygon_t> polygons) {

  // Create bbox vector TODO: rotation (cmd flag)
  std::vector<box_t> bboxes;
  bboxes.resize(polygons.size());
  std::transform(polygons.cbegin(), polygons.cend(), bboxes.begin(),
                 [](const polygon_t& p) -> box_t {
                   box_t b;
                   bg::envelope(p,b);
                   return b;
                 });

  // construct root
  s.nodes.push_back({s.bin});
  node_t* root = &s.nodes.back();

  for (box_t& bbox : bboxes)
    if (auto node = find_node(s, root, bbox)) // TODO else pack in other bin
      s.fits[&bbox] = split_node(s, node, bbox);

  // Create fitted box vector
  std::vector<box_t> fboxes;
  fboxes.resize(bboxes.size());
  std::transform(bboxes.cbegin(), bboxes.cend(), fboxes.begin(),
                 [&s](const box_t& bbox) -> box_t {
                   // TODO filter first on fit
                   node_t* n = s.fits[&bbox];
                   //if (!n) return box_t();
                   point_t min = bbox.min_corner();
                   point_t max = bbox.max_corner();
                   bg::add_point(min, n->box.min_corner());
                   bg::add_point(max, n->box.min_corner());
                   return box_t(min, max);
                 });

  return fboxes;
}

node_t* bbpack::find_node(state& s, node_t* root, box_t& bbox) {

  if (root->used) {
    node_t* up = root->up;
    node_t* right = root->right;
    if (auto node = find_node(s, right, bbox))
      return node;
    return find_node(s, up, bbox);
  }

  auto rtdims = dims(root->box);
  auto bbdims = dims(bbox);

  if ((bbdims.w <= rtdims.w) && (bbdims.h <= rtdims.h))
    return root;

  return nullptr;
}

node_t* bbpack::split_node(state& s, node_t* n, box_t& bbox) {

  s.fits[&bbox] = n;
  n->used = true;

  double n_min_x = bg::get<bg::min_corner, 0>(n->box);
  double n_min_y = bg::get<bg::min_corner, 1>(n->box);
  double n_max_x = bg::get<bg::max_corner, 0>(n->box);
  double n_max_y = bg::get<bg::max_corner, 1>(n->box);

  auto bbdims = dims(bbox);

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
