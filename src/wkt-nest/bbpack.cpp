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
  polygon_t transform(const polygon_t& p, const matrix_t& t) {
    polygon_t result = p;
    bg::for_each_point(result, [&t](point_t& p) {
                                 vector_t vec({p.x(),p.y(),1});
                                 vector_t r = ublas::prod(t, vec);
                                 bg::set<0>(p, r(0));
                                 bg::set<1>(p, r(1));
                               });
    return result;
  }
  matrix_t translation(double x, double y) {
    array_t a ({1,0,x,
                0,1,y,
                0,0,1});
    return matrix_t(3,3,std::move(a));
  }
  matrix_t translation(point_t p) { return translation(p.x(), p.y()); }
}


bbpack::item_t::item_t(const polygon_t* p) {
  _source = p;
  bg::envelope(*_source,_bbox);
  // initial transform = move to origin TODO extend with optimal rotation
  _init_transform = translation(-_bbox.min_corner().x(), -_bbox.min_corner().y());
  absolute_transform(identity_matrix(3));
}

void item_t::relative_transform(const matrix_t& t) {
  _transform = ublas::prod(_transform,t);

  // Apply transformation to source polygon
  _polygon = ::transform(*_source, _transform);
  bg::envelope(_polygon,_bbox);
}

void item_t::absolute_transform(const matrix_t& t) {
  // reset to initial transform
  _transform = _init_transform;
  relative_transform(t);
}


std::vector<matrix_t> bbpack::fit(state_t& s, const std::vector<polygon_t>& polygons) {

  for (const polygon_t& p : polygons)
    s.items.push_back(item_t(&p));

  if (s.sorting == SORTING::HEIGHT)
    std::sort(s.items.begin(), s.items.end(), [](const item_t& i1, const item_t& i2){
                                                return i1.bbox()->max_corner().y() > i2.bbox()->max_corner().y();
                                              });

  // construct root
  s.nodes.push_back({s.bin});
  node_t* root = &s.nodes.back();

  for (item_t& item : s.items)
    if (auto node = find_node(s, root, &item)) // TODO else pack in other bin
      if (split_node(s, node, &item))
        s.fits[item.source()] = &item;

  // Create the transformations vector
  std::vector<matrix_t> transformations;
  transformations.resize(s.items.size());
  std::transform(polygons.cbegin(), polygons.cend(), transformations.begin(),
                 [&s](const polygon_t& p) -> matrix_t {
                   return *(s.fits[&p]->transform());
                 });

  return transformations;
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

  n->used = true;

  double n_min_x = bg::get<bg::min_corner, 0>(n->box);
  double n_min_y = bg::get<bg::min_corner, 1>(n->box);
  double n_max_x = bg::get<bg::max_corner, 0>(n->box);
  double n_max_y = bg::get<bg::max_corner, 1>(n->box);

  item->absolute_transform(translation(n_min_x, n_min_y));

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
