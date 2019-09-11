#include "wkt-nest/wktio.hpp"

#include <sstream>
#include <boost/geometry/io/wkt/read.hpp>

using namespace wktnest;

namespace {
  const char* BOX = "BOX(";
  const char* POLYGON = "POLYGON(";

  std::optional<std::string> read_token(const char* token, std::istream &in) {
    char c;
    size_t i=0;
    // Read until token matched
    while (token[i]!=0 && in.get(c)) {
      if (c!=token[i++])
        // reset i
        i=0;
    }

    std::stringstream ps(token);
    ps.seekp(0, std::ios_base::end);

    size_t brack_level=1; // first bracket is part of the match
    while (brack_level>0 && in.get(c) && ps.put(c))
      switch (c) {
      case '(':
        brack_level++;
        break;
      case ')':
        brack_level--;
        break;
      }

    if (brack_level>0)
      return std::optional<std::string>();

    std::cerr << ps.str() << std::endl;
    return ps.str();
  }
}

namespace {
  const crd_t f=100000;

  typedef double N;
  typedef boost::geometry::model::d2::point_xy<N> point;
  typedef boost::geometry::model::polygon<point> polygon;
  typedef boost::geometry::model::box<point> box;

  template<typename C>
  C to_crd_t(const N& n) {
    return (C)(n*f);
  }
  template<>
  N to_crd_t<N>(const N& n) {
    return n*f;
  }
  inline
  crd_t to_crd_t(const N& n) {
    return to_crd_t<crd_t>(n);
  }
}

#include <boost/geometry/algorithms/for_each.hpp>
std::optional<box_t> wktnest::read_box(std::istream &in) {
  if (auto s = read_token(BOX, in)) {
    box b;
    boost::geometry::read_wkt(*s, b);
    box_t b2 = {{to_crd_t(b.min_corner().x()), to_crd_t(b.min_corner().y())},
                {to_crd_t(b.max_corner().x()), to_crd_t(b.max_corner().y())}};
    return b2;
  }
  return std::optional<box_t>();
}

std::optional<polygon_t> wktnest::read_polygon(std::istream &in) {
  if (auto s = read_token(POLYGON, in)) {
    polygon p;
    boost::geometry::read_wkt(*s, p);
    boost::geometry::for_each_point(p, [](point& pnt){
                                         pnt.x(to_crd_t(pnt.x()));
                                         pnt.y(to_crd_t(pnt.y()));
                                       });
    polygon_t p2;
    boost::geometry::transform(p, p2);
    return p2;
  }
  return std::optional<polygon_t>();
}
