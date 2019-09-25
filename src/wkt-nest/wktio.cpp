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
  if (auto s = read_token(POLYGON, in)) {
    polygon_t p;
    boost::geometry::read_wkt(*s, p);
    return p;
  }
  return std::optional<polygon_t>();
}
