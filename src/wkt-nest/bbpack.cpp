#include "wkt-nest/bbpack.hpp"

#include <algorithm>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/geometry/algorithms/for_each.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/algorithms/area.hpp>
namespace bg = boost::geometry;

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
namespace bm = boost::math;
typedef bm::policies::policy<bm::policies::evaluation_error<bm::policies::ignore_error>> ignore_eval_err;


using namespace wktnest::bbpack;
using wktnest::matrix_t;
using wktnest::identity_matrix;
using wktnest::SORTING;

typedef wktnest::N flt_t;

#define S 100000
typedef int_fast32_t crd_t;
//typedef flt_t crd_t;
const matrix_t to_bbpack = {S,0,0,
                            0,S,0,
                            0,0,1};
const matrix_t from_bbpack = {1.0/S,0,0,
                              0,1.0/S,0,
                              0,0,1};

typedef boost::geometry::model::d2::point_xy<crd_t> point_t;
typedef boost::geometry::model::polygon<point_t> polygon_t;
typedef boost::geometry::model::box<point_t> box_t;

namespace {
  const boost::uintmax_t MAX_IT = 10;
  const std::function<bool(crd_t,crd_t)> stop_condition = [](crd_t a, crd_t b) //{ return false; };
    {return std::abs(std::abs(a)-std::abs(b))<0.01;};

  struct dim_t {
    crd_t w;
    crd_t h;
  };

  dim_t dims(const box_t* b) {
    crd_t min_x = bg::get<bg::min_corner, 0>(*b);
    crd_t min_y = bg::get<bg::min_corner, 1>(*b);
    crd_t max_x = bg::get<bg::max_corner, 0>(*b);
    crd_t max_y = bg::get<bg::max_corner, 1>(*b);
    return { max_x-min_x, max_y-min_y };
  }

  matrix_t prod(const matrix_t& m1, const matrix_t& m2) {
    return m1.matrix() * m2.matrix();
  }

  polygon_t transform(const wktnest::polygon_t& p, const matrix_t& t) {
    polygon_t result;
    bg::transform(p, result, matrix_t(t.matrix() * to_bbpack.matrix()));
    return result;
  }

  matrix_t translation(double x, double y) {
    return {1,0,x,
            0,1,y,
            0,0,1};
  }
  matrix_t translation(point_t p) { return translation(p.x(), p.y()); }

  matrix_t rotation(double angle) {
    double cosa = std::cos(angle);
    double sina = std::sin(angle);
    return {cosa,sina,0,
            -sina,cosa,0,
            0,0,1};
  }

  struct node_t {
    box_t box;
    bool used = false;
    node_t* up = nullptr;
    node_t* right = nullptr;
  };

  class item_t {
    bool _placed;
    const wktnest::polygon_t* _source;
    matrix_t _init_transform;

    polygon_t _polygon;
    box_t _bbox;
    matrix_t _transform;

  public:
    item_t(const wktnest::polygon_t* p) :
      _placed(false),
      _init_transform(identity_matrix),
      _transform(identity_matrix),
      _source(p) {
      bg::envelope(*_source,_bbox);
      // initial transform = optimal rotation and move to origin
      init_transform(identity_matrix);
    }

    const wktnest::polygon_t* source() const { return _source; }
    const polygon_t* polygon() const { return &_polygon; }
    const box_t* bbox() const { return &_bbox; }
    const matrix_t* init_transform() const { return &_init_transform; }
    const matrix_t* transform() const { return &_transform; }
    void placed(bool p) { _placed = p; }
    bool placed() const { return _placed; }

    /*
     * Transformation is typically done L = T * R * S
     * first scale, then rotate, lastly
     */
    void relative_transform(const matrix_t& t) {
      _transform = prod(t, _transform);

      // Apply transformation to source polygon
      _polygon = ::transform(*_source, _transform);
      bg::envelope(_polygon,_bbox);
    }
    void absolute_transform(const matrix_t& t) {
      // reset to initial transform
      _transform = _init_transform;
      relative_transform(t);
    }
    void init_transform(const matrix_t& t) {
      absolute_transform(t);
      _init_transform = prod(translation(-_bbox.min_corner().x(), -_bbox.min_corner().y()), _transform);
      absolute_transform(identity_matrix);
    }
  };

  struct state_t {
    box_t bin;
    SORTING sorting;
    bool compact;
    std::list<item_t> items;
    std::map<const wktnest::polygon_t*, item_t*> fits;
    std::map<const item_t*, size_t> item_to_bin_idx;

    // TODO below still needed?

    // list of nodes per bin
    //std::vector<std::list<node_t>> nodes;
    // single bin of nodes TODO switch to multiple bins?
    // or better chain calls to wkt-nest?
    // Nodes should be stored in list since vector will reallocate on resizing
    std::list<node_t> nodes;
  };

  const std::vector<double> rotations = {M_PI/2};
  void create_items(state_t& s, const wktnest::polygon_t& p) {
    item_t i(&p);
    auto f_area = [&](double angle) {
                    i.absolute_transform(rotation(angle));
                    return bg::area(*i.bbox());
                  };
    boost::uintmax_t max_iter = MAX_IT;
    auto min = bm::tools::brent_find_minima(f_area, M_PI/8, M_PI*3/8, 4, max_iter);
    i.init_transform(rotation(min.first));
    s.items.push_back(i);

    for (double r : rotations) {
      item_t ir(&p);
      ir.init_transform(rotation(min.first+r));
      s.items.push_back(ir);
    }
  }
  void add_flips(std::list<item_t>& items) {
    bool before = false;
    auto it=items.begin();
    while (it!=items.end()) {
      item_t flip(it->source());
      flip.init_transform(prod(rotation(M_PI), *it->init_transform()));
      // Alternate before aft insertion (known to work well for fins)
      if (before = !before) {
        items.insert(it, flip);
        it++;
      }
      else {
        it++;
        items.insert(it, flip);
      }
    }
  }

  template<typename Geometry>
  bool within(const Geometry& g1, const Geometry& g2) {
    return bg::within(g1, g2);
  }
  template<typename Geometry>
  bool collide(const Geometry& g1, const Geometry& g2) {
    return bg::intersects(g1, g2);
  }
  template<>
  bool collide<item_t>(const item_t& i1, const item_t& i2) {
    return collide(*i1.bbox(), *i2.bbox()) && collide(*i1.polygon(), *i2.polygon());
  }
  bool can_claim_space(const item_t& item, const state_t& s) {
    // Check for collision with other items
    for (const item_t& i : s.items) {
      if (!i.placed())
        continue;
      if (collide(item, i))
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
    /* TODO make bracketing optional
    for (double d=0.0001; d<1; d=d+0.05) {
      if (std::signbit(f_collision(d))) {
        bracket.second = d;
        break;
      }
      bracket.first = d;
    }
    */
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
    std::pair<double,double> perc_move = bm::tools::bisect(f_collision, bracket.first, bracket.second, stop_condition, max_it, ignore_eval_err());

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

    auto f_perc_value = [&start_pnt](double p) -> crd_t {return p*start_pnt.x();};
    auto f_transform = [&start_pnt](crd_t v){return translation(v,start_pnt.y());};
    return find_free_space(s, item, f_perc_value, f_transform);
  }
  template<>
  bool find_free_space<DOWN>(state_t& s, item_t* item) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p) -> crd_t {return p*start_pnt.y();};
    auto f_transform = [&start_pnt](crd_t v){return translation(start_pnt.x(),v);};
    return find_free_space(s, item, f_perc_value, f_transform);
  }
  template<>
  bool find_free_space<DOWN|RIGHT>(state_t& s, item_t* item) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p) -> crd_t {return start_pnt.y()-(p*start_pnt.y());};
    auto f_transform = [&start_pnt](crd_t v){return translation(start_pnt.x()+(v/2),start_pnt.y()-v);};
    return find_free_space(s, item, f_perc_value, f_transform);
  }
}

node_t* find_node(state_t& s, node_t* root, item_t* item, size_t rec_depth=0);

fit_result wktnest::bbpack::fit(const wktnest::box_t& bin, const std::vector<wktnest::polygon_t>& polygons, SORTING sorting, bool compact) {
  ::box_t lb;
  bg::transform(bin, lb, to_bbpack);
  state_t s = { lb, sorting, compact };

  for (const wktnest::polygon_t& p : polygons)
    create_items(s, p);

  if (s.sorting == SORTING::HEIGHT)
    s.items.sort([](const item_t& i1, const item_t& i2){
                   return i1.bbox()->max_corner().y() > i2.bbox()->max_corner().y();
                 });

  add_flips(s.items);

  // construct root
  s.nodes.push_back({s.bin});
  node_t* root = &s.nodes.back();

  for (item_t& item : s.items)
    if (!s.fits[item.source()])
      if (auto node = find_node(s, root, &item)) // TODO else pack in other bin
        s.fits[item.source()] = &item;

  // Create the result vector
  std::vector<placement> result;
  result.resize(polygons.size());
  std::transform(polygons.cbegin(), polygons.cend(), result.begin(),
                 [&s](const polygon_t& p) -> placement {
                   if (!s.fits[&p])
                     return {0};
                   item_t* item = s.fits[&p];
                   wktnest::polygon_t rp;
                   wktnest::box_t rb;
                   bg::transform(*item->polygon(), rp, from_bbpack);
                   bg::transform(*item->bbox(), rb, from_bbpack);
                   return {1,rp,rb,*item->transform()};
                 });

  return result;
}

node_t* split_node(state_t& s, node_t* node, item_t* item) {

  crd_t n_min_x = bg::get<bg::min_corner, 0>(node->box);
  crd_t n_min_y = bg::get<bg::min_corner, 1>(node->box);
  crd_t n_max_x = bg::get<bg::max_corner, 0>(node->box);
  crd_t n_max_y = bg::get<bg::max_corner, 1>(node->box);

  item->absolute_transform(translation(n_min_x, n_min_y));

  if (s.compact) {
    for (size_t i=0; i<MAX_IT; i++) {
      bool down = find_free_space<DOWN|RIGHT>(s, item);
      bool left = find_free_space<LEFT>(s, item);
      if (!left && !down)
        // TODO apparently never reached?
        break;
    }

    if (!can_claim_space(*item,s))
      return nullptr;
  }

  node->used = true;
  item->placed(true);

  /*
  // if bbox outside node, free up full node
  if (!collide(*item->bbox(), node->box)) {
    // node up
    node->up = nullptr;
    // node right
    box_t right = node->box;
    s.nodes.push_back({ right });
    node->right = &s.nodes.back();
    return node;
  }
  */

  /* Splitting on bbox is more robust when having many equal sized objects
   * A small downwards move of an element can make the entire row unusable
  crd_t x_split = std::max(n_min_x, std::min(std::ceil(item->bbox()->max_corner().x()), n_max_x));
  crd_t y_split = std::max(n_min_y, std::min(std::ceil(item->bbox()->max_corner().y()), n_max_y));
  */
  crd_t x_split = std::max(n_min_x, std::min(item->bbox()->max_corner().x(), n_max_x));
  crd_t y_split = std::max(n_min_y, std::min(item->bbox()->max_corner().y(), n_max_y));

  /*
  crd_t x_split = std::ceil(n_min_x + dims(item->bbox()).w);
  crd_t y_split = std::ceil(n_min_y + dims(item->bbox()).h);

  auto item_box = item->bbox();
  crd_t x_split = item_box->max_corner().x();
  crd_t y_split = item_box->max_corner().y();
  */

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

node_t* find_node(state_t& s, node_t* root, item_t* item, size_t rec_depth) {

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

  auto fit_factor = s.compact? 2 : 1;
  // do fit height with fit_factor as compaction might still push it inside
  // do not use fit_factor for width, since non-fits will be skipped, while there might be place up
  if ((bbdims.w <= rtdims.w) && (bbdims.h <= rtdims.h*fit_factor))
    return split_node(s, root, item);

  return nullptr;
}
