#include "wkt-nest/bbpack-geometry.hpp"

#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/centroid.hpp>
#include <boost/geometry/strategies/cartesian/centroid_bashein_detmer.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/algorithms/area.hpp>
namespace bg = boost::geometry;

using namespace bbpack::geometry;

flt_t bbpack::geometry::area(const box_t& b) {
  return bg::area(b);
}

void bbpack::geometry::centroid(const box_t& b, point_t& c) {
  return bg::centroid(b, c);
}
const bg::strategy::centroid::bashein_detmer<point_t,point_t,crd_t> centroid_strategy;
void bbpack::geometry::centroid(const polygon_t& b, point_t& c) {
  return bg::centroid(b, c, centroid_strategy);
}

bool bbpack::geometry::intersects(const box_t& i1, const box_t& i2) {
  return bg::intersects(i1, i2);
}
bool bbpack::geometry::intersects(const box_t& i1, const polygon_t& i2) {
  return bg::intersects(i1, i2);
}
bool bbpack::geometry::intersects(const polygon_t& i1, const polygon_t& i2) {
  return bg::intersects(i1, i2);
}
