#include "wkt-nest/wktio.hpp"

#include <sstream>
#include <boost/geometry/io/wkt/read.hpp>
#include <wkt-nest/bbpack-geometry.hpp>

using namespace wktnest;

namespace {
  const char* BOX = "BOX(";
  const char* POLYGON = "POLYGON(";
  const char* MULTIPOLYGON = "MULTIPOLYGON(";

  std::optional<std::string> read_token(const char* token, std::istream &in) {
    char c;
    size_t i=0;
    // Read until token matched
    // Or bail out if char differs
    while (token[i]!=0) {
      c=in.peek();
      if (!std::isspace(c) && c!=token[i++])
        return std::optional<std::string>();
      in.get(c);
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

    return ps.str();
  }
}

std::optional<box_t> wktnest::read_box(std::istream &in) {
  if (auto s = read_token(BOX, in)) {
    box_t b;
    boost::geometry::read_wkt(*s, b);
    return b;
  }
  return std::optional<box_t>();
}

std::optional<polygon_t> wktnest::read_polygon(std::istream &in) {
  if (auto s = read_token(MULTIPOLYGON, in)) {
    multi_polygon_t mp;
    boost::geometry::read_wkt(*s, mp);
    // TODO more?
    return bbpack::geometry::union_(mp[0],mp[1]);
  }
  if (auto s = read_token(POLYGON, in)) {
    polygon_t p;
    boost::geometry::read_wkt(*s, p);
    return p;
  }
  return std::optional<polygon_t>();
}
