#include "wkt-nest/bbpack.hpp"
#include "wkt-nest/bbpack-geometry.hpp"

#include <iostream>
#include <algorithm>
#include <boost/math/special_functions/pow.hpp>
#include <boost/geometry/algorithms/envelope.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/algorithms/transform.hpp>
namespace bg = boost::geometry;

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
namespace bm = boost::math;
typedef bm::policies::policy<bm::policies::evaluation_error<bm::policies::ignore_error>> ignore_eval_err;

using namespace wktnest::bbpack;
using namespace bbpack::geometry;
using wktnest::SORTING;

namespace {
  const boost::uintmax_t MAX_IT = 10;
  const std::function<bool(flt_t,flt_t)> stop_condition = [](flt_t a, flt_t b) //{ return false; };
    {return std::abs(std::abs(a)-std::abs(b))<0.01;};

  template<typename T>
  struct S
  {
    static constexpr flt_t value = 65536;// bm::pow<16>(2);
  };
  template<>
  struct S<flt_t>
  {
    static constexpr flt_t value = 1;
  };

  const matrix_t m_to_bbpack = {S<crd_t>::value,0,0,
                              0,S<crd_t>::value,0,
                              0,0,1};

  template<typename T>
  typename std::enable_if<std::negation<std::is_same<T, polygon_t>>::value, polygon_t>::type
  to_bbpack(const T& p, const matrix_t& t) {
    polygon_t result;
    bg::transform(p, result, matrix_t(t.matrix() * m_to_bbpack.matrix()));
    return result;
  }
  template<typename T>
  typename std::enable_if<std::negation<std::is_same<T, box_t>>::value, box_t>::type
  to_bbpack(const T& b) {
    box_t result;
    bg::transform(b, result, m_to_bbpack);
    return result;
  }
  polygon_t to_bbpack(const polygon_t& p, const matrix_t& t) {
    polygon_t result;
    bg::transform(p, result, t);
    return result;
  }
  box_t to_bbpack(const box_t& b) {
    return b;
  }

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
    const wktnest::polygon_t* _source;
    polygon_t _buffered;
    matrix_t _init_transform;

    polygon_t _polygon;
    box_t _bbox;
    matrix_t _transform;

    flt_t _footprint;
    point_t _centroid;

  public:
    int rand;

    // Only allow moving, no copying
    item_t(item_t&&) = default;
    item_t() = delete;
    void operator=(item_t const &x) = delete;

    item_t(const wktnest::polygon_t* p, flt_t buffer_distance) :
      _init_transform(identity_matrix),
      _transform(identity_matrix),
      _source(p)
    {
      _buffered = buffer(to_bbpack(*_source, _transform), buffer_distance);
      // initial transform = optimal rotation and move to origin
      init_transform(identity_matrix);
      // random number for shuffling
      rand = std::rand();
      _footprint = area(_buffered);
    }

  private:
    item_t(const item_t& other) :
      _init_transform(identity_matrix),
      _transform(identity_matrix),
      _source(other._source), _buffered(other._buffered),
      _footprint(other._footprint)
      {
        // initial transform = optimal rotation and move to origin
        init_transform(identity_matrix);
        // random number for shuffling
        rand = std::rand();
      }
  public:
    item_t clone() { return item_t(*this); }

    const wktnest::polygon_t* source() const { return _source; }
    const polygon_t* polygon() const { return &_polygon; }
    const box_t* bbox() const { return &_bbox; }
    const matrix_t* init_transform() const { return &_init_transform; }
    const matrix_t* transform() const { return &_transform; }
    flt_t footprint() const { return _footprint; }
    point_t centroid() const { return _centroid; }

    /*
     * Transformation is typically done L = T * R * S
     * first scale, then rotate, lastly
     */
    void relative_transform(const matrix_t& t) {
      _transform = prod(t, _transform);

      // Apply transformation to buffered polygon
      _polygon = to_bbpack(_buffered, _transform);
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
      ::centroid(_polygon, _centroid);
    }
  };

  struct state_t {
    box_t bin;
    flt_t buffer_distance;
    SORTING sorting;
    bool compact;
    bool centroid_left = true;

    // or better chain calls to wkt-nest?
    // Nodes should be stored in list since vector will reallocate on resizing
    std::list<node_t> nodes;

    //const std::vector<double> rotations = {M_PI, M_PI/2, M_PI/4, M_PI*3/4, M_PI/2+M_PI, M_PI/4+M_PI, M_PI*3/4+M_PI};
    const std::vector<double> rotations = {M_PI, M_PI/2, M_PI/2+M_PI};
    //const std::vector<double> rotations = {M_PI};

  private:
    box_t _union_box = {{0,0},{0,0}};
    polygon_t _union_poly;
    std::list<item_t> _item_cont;
    std::list<item_t*> _source_items;
    std::map<const wktnest::polygon_t*, std::vector<item_t*>> items_for;
    std::map<const item_t*, size_t> item_to_bin_idx;
  public:
    state_t (box_t bin, flt_t buffer_distance, SORTING sorting, bool compact)
      : bin(bin), buffer_distance(buffer_distance), sorting(sorting), compact(compact) {}

    const std::list<item_t>& all_items() const { return _item_cont; }
    std::list<item_t*>& items() { return _source_items; }
    const std::list<item_t*>& items() const { return _source_items; }
    void add_item(item_t&& item) {
      _item_cont.push_back(std::move(item));
      item_t* i = &_item_cont.back();
      _source_items.push_back(i);

      std::vector<item_t*>& mutations = items_for[i->source()];
      mutations.push_back(i);
      item_to_bin_idx[i] = 0;

      matrix_t itm = i->transform()->matrix();
      for (double r : rotations) {
        item_t ir(i->source(), buffer_distance);
        ir.init_transform(matrix_t(itm.matrix() * rotation(r).matrix()));
        _item_cont.push_back(std::move(ir));
        item_to_bin_idx[&_item_cont.back()] = 0;
        mutations.push_back(&_item_cont.back());
      }
    }
    void place(const item_t& item) {
      item_to_bin_idx[&item] = 1;
      _union_box = union_(_union_box, *item.bbox());
      _union_poly = union_(_union_poly, convert(*item.bbox()));
    };
    const std::vector<item_t*>& mutations(const wktnest::polygon_t* source) const {
      return items_for.at(source);
    }
    const std::vector<item_t*>& mutations(const item_t& item) const {
      return items_for.at(item.source());
    }
    bool placed (const item_t& item) const {
      return item_to_bin_idx.at(&item) > 0;
    }
    item_t* placed (const wktnest::polygon_t* source) const {
      const std::vector<item_t*>& its = mutations(source);
      auto it = std::find_if(its.cbegin(), its.cend(), [this](const item_t* i){ return placed(*i); });
      if (it == its.end())
        return nullptr;
      return *it;
    }

    box_t& union_box() { return _union_box; }
    polygon_t& union_poly() { return _union_poly; }
  };

  void create_items(state_t& s, const wktnest::polygon_t& p) {
    item_t i(&p, s.buffer_distance);
    auto f_area = [&](double angle) {
                    i.absolute_transform(rotation(angle));
                    return dims(i.bbox()).w; // minimal width = maximal height
                  };
    boost::uintmax_t max_iter = MAX_IT;
    auto min = bm::tools::brent_find_minima(f_area, M_PI/8, M_PI*7/8, 4, max_iter);
    i.init_transform(rotation(min.first));
    s.add_item(std::move(i));
  }

  template<typename Geometry1, typename Geometry2>
  bool within(const Geometry1& g1, const Geometry2& g2) {
    return bg::within(g1, g2);
  }
  template<>
  bool within<item_t,box_t>(const item_t& g1, const box_t& g2) {
    return bg::within(*g1.bbox(), g2);
  }
  template<typename Geometry1, typename Geometry2>
  bool collide(const Geometry1& g1, const Geometry2& g2) {
    return intersects(g1, g2);
  }
  template<>
  bool collide<box_t, box_t>(const box_t& i1, const box_t& i2) {
    // Override box collision to ensure touching does not count as collision
    point_t ca, cb;
    centroid(i1, ca);
    centroid(i2, cb);
    auto da = dims(&i1);
    auto db = dims(&i2);
    return (std::abs(ca.x() - cb.x()) * 2 < (da.w + db.w)) &&
           (std::abs(ca.y() - cb.y()) * 2 < (da.h + db.h));
  }
  template<>
  bool collide<item_t, item_t>(const item_t& i1, const item_t& i2) {
    return collide(*i1.bbox(), *i2.bbox()) && collide(*i1.polygon(), *i2.polygon());
  }
  template<>
  bool collide<box_t, item_t>(const box_t& i1, const item_t& i2) {
    return collide(i1, *i2.bbox()) && collide(i1, *i2.polygon());
  }

  template<typename Geometry>
  bool can_claim_space(const Geometry& item, const state_t& s) {
    // Check for collision with other items
    // TODO get placed items from state?
    for (const item_t& i : s.all_items()) {
      if (!s.placed(i))
        continue;
      if (collide(item, i))
        return false;
    }
    return within(item, s.bin);
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
  std::pair<double, double> bracket_solution(const std::function<double(double)>& f_collision,
                                             std::pair<double,double> bracket = {0,1}) {
    double min = bracket.first;
    double max = bracket.second;
    for (double d=min; d<max; d=d+0.01) {
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
                       const std::function<matrix_t(double)>& f_transform,
                       bool bracket=false) {
    matrix_t initial_transform = *item->transform();
    // Move to left until collision
    auto f_collision = [&](double perc) -> double {
                         if (perc==0) return 1;
                         double v = f_perc_value(perc);
                         item->absolute_transform(f_transform(v));
                         return can_claim_space(*item, s)? -1*perc : 1*perc;
                       };

    std::pair<double, double> range = bracket? bracket_solution(f_collision) : std::pair<double, double>(0,1);

    boost::uintmax_t max_it = MAX_IT;
    std::pair<double,double> perc_move = bm::tools::bisect(f_collision, range.first, range.second, stop_condition, max_it, ignore_eval_err());

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
  bool find_free_space(state_t& s, item_t* item, bool bracket=false) {
    assert(false);
    return false;
  }
  template<>
  bool find_free_space<LEFT>(state_t& s, item_t* item, bool bracket) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p) -> crd_t {return p*start_pnt.x();};
    auto f_transform = [&start_pnt](crd_t v){return translation(v,start_pnt.y());};
    return find_free_space(s, item, f_perc_value, f_transform, bracket);
  }
  template<>
  bool find_free_space<DOWN>(state_t& s, item_t* item, bool bracket) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p) -> crd_t {return p*start_pnt.y();};
    auto f_transform = [&start_pnt](crd_t v){return translation(start_pnt.x(),v);};
    return find_free_space(s, item, f_perc_value, f_transform, bracket);
  }
  template<>
  bool find_free_space<DOWN|RIGHT>(state_t& s, item_t* item, bool bracket) {
    point_t start_pnt = item->bbox()->min_corner();

    auto f_perc_value = [&start_pnt](double p) -> crd_t {return start_pnt.y()-(p*start_pnt.y());};
    auto f_transform = [&start_pnt](crd_t v){return translation(start_pnt.x()+(v*2),start_pnt.y()-v);};
    return find_free_space(s, item, f_perc_value, f_transform, bracket);
  }

  inline
  bool can_grow_up(node_t* node) {
    if (node->up && node->up->used)
      return false;
    return node->box.min_corner().x() == 0;
  }
}

node_t* find_node(state_t& s, node_t* root, item_t* item, size_t rec_depth=0);

fit_result wktnest::bbpack::fit(const wktnest::box_t& bin, const std::vector<wktnest::polygon_t>& polygons, double distance, SORTING sorting, bool compact) {
  std::srand(std::time(0));

  state_t s = { to_bbpack(bin), distance*S<crd_t>::value/2, sorting, compact };

  for (const wktnest::polygon_t& p : polygons)
    create_items(s, p);

  switch (s.sorting) {
  case SORTING::HEIGHT:
    s.items().sort([](const item_t* i1, const item_t* i2){
                   return i1->bbox()->max_corner().y() > i2->bbox()->max_corner().y();
                 });
    break;
  case SORTING::AREA:
    s.items().sort([](const item_t* i1, const item_t* i2){
                   if (i1->footprint() == i2->footprint())
                     return i1->bbox()->max_corner().y() > i2->bbox()->max_corner().y();
                   return i1->footprint() > i2->footprint();
                 });
    break;
  case SORTING::SHUFFLE:
    s.items().sort([](const item_t* i1, const item_t* i2){
                   return i1->rand < i2->rand;
                 });
  }

  // construct root
  auto id = dims(s.items().front()->bbox());
  s.nodes.push_back({ {{0,0},{s.bin.max_corner().x(),id.h}} });
  node_t* root = &s.nodes.back();

  for (item_t* item : s.items())
    if (!s.placed(item->source()))
      find_node(s, root, item);

  // Create the result vector
  std::vector<placement> result;
  result.resize(polygons.size());
  std::transform(polygons.cbegin(), polygons.cend(), result.begin(),
                 [&s](const polygon_t& p) -> placement {
                   if (!s.placed(&p))
                     return {0};
                   item_t* item = s.placed(&p);
                   auto t = item->transform()->matrix();
                   t.a[0][2] /= S<crd_t>::value;
                   t.a[1][2] /= S<crd_t>::value;
                   wktnest::polygon_t rp;
                   bg::transform(*item->source(), rp, matrix_t(t));
                   return {1,
                           rp,
                           t};
                 });

  // Output the nesting efficiency
  flt_t total_footprint = 0;
  for (const polygon_t& p : polygons)
    if (item_t* item = s.placed(&p))
      total_footprint += item->footprint();
  std::cerr << "union footprint: " << total_footprint/area(s.union_box()) * 100 << "%" << std::endl;
  std::cerr << "box footprint:   " << total_footprint/(dims(&s.bin).w * dims(&s.union_box()).h) * 100 << "%" << std::endl;
  std::cerr << "bin footprint:   " << total_footprint/area(s.bin) * 100 << "%" << std::endl;

  return result;
}

bool split_node(state_t& s, node_t* node, item_t* item, bool artificial = false) {

  crd_t n_min_x = bg::get<bg::min_corner, 0>(node->box);
  crd_t n_min_y = bg::get<bg::min_corner, 1>(node->box);
  crd_t n_max_x = bg::get<bg::max_corner, 0>(node->box);
  crd_t n_max_y = bg::get<bg::max_corner, 1>(node->box);

  item->absolute_transform(translation(n_min_x, n_min_y));

  if (s.compact) {
    // If new row, try placing in top right first (might still fit row below)
    if (n_min_x==0 && n_min_y!=0 && !artificial) {
      auto id = dims(item->bbox());
      // Artificial top right node
      node_t artificial_node = { {{n_max_x-id.w,n_max_y-id.h},{n_max_x, n_max_y}} };
      if (split_node(s, &artificial_node, item, true))
        return true;
    }

    auto mutations = s.mutations(*item);
    std::sort(mutations.begin(), mutations.end(), [&s](const item_t* i1, const item_t* i2){
                                                    if (dims(i1->bbox()).h > dims(i2->bbox()).h*1.1)
                                                      return true;
                                                    if (dims(i2->bbox()).h > dims(i1->bbox()).h*1.1)
                                                      return false;
                                                    // alternating centroid sorting
                                                    bool right = i1->centroid().x() > i2->centroid().x();
                                                    return s.centroid_left != right;
                                                  });
    crd_t prev_dist = std::numeric_limits<crd_t>::max();
    flt_t prev_area = std::numeric_limits<flt_t>::max();
    flt_t ref_area = area(s.union_poly());
    item = mutations.front();
    int max_mut = 2;
    for (item_t* mut : mutations) {
      mut->absolute_transform(translation(n_min_x, n_min_y));
      if (mut->bbox()->max_corner().x() > n_max_x)
        mut->relative_transform(translation(n_max_x-mut->bbox()->max_corner().x(),0));
      // Find free space
      // First push left & down using bracketing
      find_free_space<LEFT>(s, mut, true);
      find_free_space<DOWN>(s, mut, true);
      // Then iterate zigzag (LEFT & DOWN|RIGHT)
      for (size_t i=0; i<MAX_IT; i++) {
        bool left = find_free_space<LEFT>(s, mut);
        bool down = find_free_space<DOWN|RIGHT>(s, mut);
        if (!left && !down)
          break;
      }
      if (!can_claim_space(*mut,s))
        // Bail out when we failed to find free space
        continue;

      crd_t new_dist = bm::pow<2>(mut->bbox()->min_corner().x()) + bm::pow<2>(mut->bbox()->min_corner().y());
      flt_t new_area = union_area(s.union_poly(), *mut->polygon())-ref_area;
      //if (new_dist<prev_dist*0.75) {
      if (new_area<prev_area*0.75) {
        prev_dist = new_dist;
        prev_area = new_area;
        item = mut;
      }

      if (--max_mut < 0)
        // Do not evaluate more fitting mutations than max_mut
        break;
    }

    if (!can_claim_space(*item,s))
      // Bail out when we failed to find free space
      return false;
  }

  node->used = true;
  s.place(*item);
  s.centroid_left = !s.centroid_left;

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

  /*
  crd_t x_split = std::max(n_min_x, std::min(std::ceil(item->bbox()->max_corner().x()), n_max_x));
  crd_t y_split = std::max(n_min_y, std::min(std::ceil(item->bbox()->max_corner().y()), n_max_y));
  */
  crd_t x_split = std::max(n_min_x, std::min(item->bbox()->max_corner().x(), n_max_x));
  crd_t y_split = std::max(n_min_y, std::min(item->bbox()->max_corner().y(), n_max_y));
  /*
  crd_t x_split = std::ceil(n_min_x + dims(item->bbox()).w);
  crd_t y_split = std::ceil(n_min_y + dims(item->bbox()).h);
  */

  // node up
  if (!can_grow_up(node)) {
    // leave creating the up node to the grow_up function
    box_t up = {{n_min_x, y_split}, {n_max_x, n_max_y}};
    s.nodes.push_back({ up });
    node->up = &s.nodes.back();
  }

  // node right
  box_t right = {{x_split, n_min_y}, {n_max_x, y_split}};
  s.nodes.push_back({ right });
  node->right = &s.nodes.back();

  return true;
}

node_t* grow_up(state_t& s, node_t* root, item_t* item, size_t rec_depth) {
  auto bin_limit = s.bin.max_corner();

  // Do not create a node on top of the top node
  if (root->box.max_corner().y() >= bin_limit.y())
    return nullptr;

  crd_t item_height = dims(item->bbox()).h;
  crd_t bot = s.union_box().max_corner().y();
  crd_t top = std::min(bot + item_height, bin_limit.y());

  if (s.compact && top >= bin_limit.y()) {
    // The top node should be fully inside the bin
    // to ensure that a valid search is performed
    top = bin_limit.y();
    bot = top - item_height;
  }

  s.nodes.push_back({ {{0, bot},
                       {bin_limit.x(), top} }});
  node_t* up = &s.nodes.back();
  root->up = up;

  s.centroid_left = true;
  if (node_t* node = find_node(s, up, item, ++rec_depth))
    return node;
  else
    return nullptr;
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
    if (can_grow_up(root))
      if (node_t* node = grow_up(s, root, item, rec_depth))
        return node;
    return nullptr;
  }

  auto rtdims = dims(&root->box);
  auto bbdims = dims(item->bbox());

  auto fit_factor = s.compact? 2 : 1;
  // do fit height with fit_factor as compaction might still push it inside
  // do not use fit_factor for width, since non-fits will be skipped, while there might be place up
  if ((bbdims.w <= rtdims.w*fit_factor) && (bbdims.h <= rtdims.h*fit_factor))
    if (split_node(s, root, item))
      return root;

  return nullptr;
}
