#include "wkt-nest/bbpack.hpp"

#include <algorithm>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/geometry/algorithms/for_each.hpp>
#include <boost/geometry/algorithms/overlaps.hpp>
#include <boost/geometry/algorithms/within.hpp>
namespace bg = boost::geometry;

#include <boost/math/tools/roots.hpp>
namespace roots = boost::math::tools;
typedef boost::math::policies::policy<boost::math::policies::evaluation_error<boost::math::policies::ignore_error>> ignore_eval_err;


using namespace wktnest;
using namespace bbpack;

namespace {
  const boost::uintmax_t MAX_IT = 10;
  const std::function<bool(double,double)> stop_condition = [](double a, double b) //{ return false; };
    {return std::abs(std::abs(a)-std::abs(b))<0.01;};

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
  matrix_t rotation(double angle) {
    matrix_t rot = identity_matrix(3);
    rot(0,0) = std::cos(angle);
    rot(1,0) = std::sin(angle);
    rot(0,1) = -rot(1,0);
    rot(1,1) = rot(0,0);
    return rot;
  }
}


bbpack::item_t::item_t(const polygon_t* p) :
  _placed(false),
  _init_transform(identity_matrix(3)),
  _transform(identity_matrix(3)),
  _source(p) {
  bg::envelope(*_source,_bbox);
  // initial transform = optimal rotation and move to origin
  init_transform(identity_matrix(3));
}

void item_t::init_transform(const matrix_t& t) {
  absolute_transform(t);
  _init_transform = ublas::prod(translation(-_bbox.min_corner().x(), -_bbox.min_corner().y()), _transform);
  absolute_transform(identity_matrix(3));
}

void item_t::relative_transform(const matrix_t& t) {
  _transform = ublas::prod(t, _transform);

  // Apply transformation to source polygon
  _polygon = ::transform(*_source, _transform);
  bg::envelope(_polygon,_bbox);
}

void item_t::absolute_transform(const matrix_t& t) {
  // reset to initial transform
  _transform = _init_transform;
  relative_transform(t);
}


#include <boost/math/tools/minima.hpp>
#include <boost/geometry/algorithms/area.hpp>
namespace {
  const std::vector<double> rotations = {M_PI/2, M_PI, M_PI*3/2};
  void create_items(state_t& s, const polygon_t& p) {
    item_t i(&p);
    auto f_area = [&](double angle) -> double {
                    i.absolute_transform(rotation(angle));
                    return bg::area(*i.bbox());
                  };
    boost::uintmax_t max_iter = MAX_IT;
    auto min = boost::math::tools::brent_find_minima(f_area, M_PI/8, M_PI*3/8, 4, max_iter);
    i.init_transform(rotation(min.first));
    s.items.push_back(i);

    for (double r : rotations) {
      item_t ir(&p);
      ir.init_transform(rotation(min.first+r));
      s.items.push_back(ir);
    }
  }
}

node_t* bbpack::find_node(state_t& s, node_t* root, item_t* item, size_t rec_depth) {

  if (rec_depth>255)
    return nullptr;
  if (root->used) {
    node_t* up = root->up;
    node_t* right = root->right;
    if (right)
      if (node_t* node = find_node(s, right, item, ++rec_depth))
        return node;
    if (up)
      if (node_t* node = find_node(s, up, item, ++rec_depth))
        return node;
    return nullptr;
  }

  auto rtdims = dims(&root->box);
  auto bbdims = dims(item->bbox());
  if ((bbdims.w <= rtdims.w*2) && (bbdims.h <= rtdims.h*2))
    return root;

  return nullptr;
}

namespace {
  template<typename Geometry>
  bool within(const Geometry& g1, const Geometry& g2) {
    return bg::within(g1, g2);
  }
  template<typename Geometry>
  bool overlaps(const Geometry& g1, const Geometry& g2) {
    return bg::overlaps(g1, g2) || bg::within(g1, g2);
  }
  template<>
  bool overlaps<item_t>(const item_t& i1, const item_t& i2) {
    return overlaps(*i1.bbox(), *i2.bbox()) && overlaps(*i1.polygon(), *i2.polygon());
  }
  bool can_claim_space(const item_t& item, const state_t& s) {
    // Check for collision with other items
    for (const item_t& i : s.items) {
      if (!i.placed())
        continue;
      if (overlaps(item, i))
        return false;
    }
    return within(*item.bbox(), s.bin);
  }

  enum direction {
                  ANY   =0,
                  LEFT  =0x1,
                  RIGHT =0x2,
                  UP    =0x4,
                  DOWN  =0x8,

                  FALLING = LEFT | DOWN,
                  RISING  = RIGHT | UP,
  };
  std::pair<double, double> bracket_solution(const std::function<double(double)>& f_collision) {
    std::pair<double,double> bracket = {0,1};
    for (double d=0.0001; d<1; d=d+0.05) {
      if (std::signbit(f_collision(d))) {
        bracket.second = d;
        break;
      }
      bracket.first = d;
    }
    return bracket;
  }
  bool find_free_space(state_t& s, item_t* item,
                       const std::function<double(double)>& f_perc_value,
                       const std::function<matrix_t(double)>& f_transform) {
    matrix_t initial_transform = *item->transform();
    // Move to left until collision
    auto f_collision = [&](double perc) -> double {
                         if (perc==0) return 1;
                         double v = f_perc_value(perc);
                         item->absolute_transform(f_transform(v));
                         return can_claim_space(*item, s)? -1*perc : 1*perc;
                       };

    std::pair<double, double> bracket = bracket_solution(f_collision);

    boost::uintmax_t max_it = MAX_IT;
    std::pair<double,double> perc_move = roots::bisect(f_collision, bracket.first, bracket.second, stop_condition, max_it, ignore_eval_err());

    if (perc_move.first<0 || perc_move.second<0) {
      // reset to start point
      item->absolute_transform(initial_transform);
      return false;
    }

    double free_perc = f_collision(perc_move.first)<0? perc_move.first : perc_move.second;
    item->absolute_transform(f_transform(f_perc_value(free_perc)));
    return free_perc<1;
  }
  template<size_t DIRECTION>
  bool find_free_space(state_t& s, item_t* item) {
    assert(false);
    return false;
  }
  template<>
  bool find_free_space<LEFT>(state_t& s, item_t* item) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p){return p*start_pnt.x();};
    auto f_transform = [&start_pnt](double v){return translation(v,start_pnt.y());};
    return find_free_space(s, item, f_perc_value, f_transform);
  }
  template<>
  bool find_free_space<DOWN>(state_t& s, item_t* item) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p){return p*start_pnt.y();};
    auto f_transform = [&start_pnt](double v){return translation(start_pnt.x(),v);};
    return find_free_space(s, item, f_perc_value, f_transform);
  }
}

std::vector<matrix_t> bbpack::fit(state_t& s, const std::vector<polygon_t>& polygons) {

  s.items.reserve(polygons.size());
  for (const polygon_t& p : polygons)
    create_items(s, p);

  if (s.sorting == SORTING::HEIGHT)
    std::sort(s.items.begin(), s.items.end(), [](const item_t& i1, const item_t& i2){
                                                return i1.bbox()->max_corner().y() > i2.bbox()->max_corner().y();
                                              });

  // construct root
  s.nodes.push_back({s.bin});
  node_t* root = &s.nodes.back();

  for (item_t& item : s.items)
    if (!s.fits[item.source()])
      if (auto node = find_node(s, root, &item)) // TODO else pack in other bin
        if (split_node(s, node, &item))
          s.fits[item.source()] = &item;

  // Create the transformations vector
  std::vector<matrix_t> transformations;
  transformations.resize(s.items.size());
  std::transform(polygons.cbegin(), polygons.cend(), transformations.begin(),
                 [&s](const polygon_t& p) -> matrix_t {
                   if (!s.fits[&p])
                     return identity_matrix(3);
                   return *(s.fits[&p]->transform());
                 });

  return transformations;
}

node_t* bbpack::split_node(state_t& s, node_t* node, item_t* item) {

  double n_min_x = bg::get<bg::min_corner, 0>(node->box);
  double n_min_y = bg::get<bg::min_corner, 1>(node->box);
  double n_max_x = bg::get<bg::max_corner, 0>(node->box);
  double n_max_y = bg::get<bg::max_corner, 1>(node->box);

  item->absolute_transform(translation(n_min_x, n_min_y));

  if (s.compact) {
    for (size_t i=0; i<MAX_IT; i++)
      if (!find_free_space<LEFT>(s, item) && !find_free_space<DOWN>(s, item))
        // TODO apparently never reached?
        break;

    if (!can_claim_space(*item,s))
      return nullptr;
  }

  node->used = true;
  item->placed(true);

  // if bbox outside node, free up full node
  if (!overlaps(*item->bbox(), node->box)) {
    // node up
    node->up = nullptr;
    // node right
    box_t right = node->box;
    s.nodes.push_back({ right });
    node->right = &s.nodes.back();
    return node;
  }

  /* Splitting on bbox is more robust when having many equal sized objects
   * A small downwards move of an element can make the entire row unusable
  double x_split = std::max(n_min_x, std::min(std::ceil(item->bbox()->max_corner().x()), n_max_x));
  double y_split = std::max(n_min_y, std::min(std::ceil(item->bbox()->max_corner().y()), n_max_y));
  /*
  */
  double x_split = std::ceil(n_min_x + dims(item->bbox()).w);
  double y_split = std::ceil(n_min_y + dims(item->bbox()).h);

  // node up
  box_t up = {{n_min_x, y_split}, {n_max_x, n_max_y}};
  s.nodes.push_back({ up });
  node->up = &s.nodes.back();

  // node right
  box_t right = {{x_split, n_min_y}, {n_max_x, y_split}};
  s.nodes.push_back({ right });
  node->right = &s.nodes.back();

  return node;
}
