#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace py = pybind11;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


class StraightSkeleton2 {
  const boost::shared_ptr<CGAL::Straight_skeleton_2<K>> _ref;

public:
  StraightSkeleton2(const boost::shared_ptr<CGAL::Straight_skeleton_2<K>> &ref) : _ref(ref) {
  }

  py::array_t<double> bisectors() const {
    size_t n = 0;
    for (auto i = _ref->halfedges_begin(); i != _ref->halfedges_end(); i++) {
      if (i->is_bisector()) {
        n += 1;
      }
    }    

    py::array_t<double> pts({ n, size_t(2), size_t(2) });
    auto r = pts.mutable_unchecked<3>();
    size_t j = 0;

    for (auto i = _ref->halfedges_begin(); i != _ref->halfedges_end(); i++) {
      if (i->is_bisector()) {
        const auto &pt1 = i->opposite()->vertex()->point();
        r(j, 0, 0) = pt1.x();
        r(j, 0, 1) = pt1.y();

        const auto &pt2 = i->vertex()->point();
        r(j, 1, 0) = pt2.x();
        r(j, 1, 1) = pt2.y();

        j += 1;
      }
    }

    assert(j == n);

    return pts;
  }
};

typedef std::shared_ptr<StraightSkeleton2> StraightSkeleton2Ref;



class Polygon2 {
  CGAL::Polygon_2<K> _poly;

public:
  Polygon2(py::array_t<double> p_vertices) {
    auto r = p_vertices.unchecked<2>();
    if (r.shape(1) != 2) {
      throw std::runtime_error("vertices need to be 2 dimensional");
    }
    const ssize_t n = r.shape(0);

    for (ssize_t i = 0; i < n; i++) {
      _poly.push_back(K::Point_2(r(i, 0), r(i, 1)));
    }
  }

  inline const CGAL::Polygon_2<K> &cgal_polygon_2() const {
    return _poly;
  }

  StraightSkeleton2Ref create_interior_straight_skeleton() const {
    return std::make_shared<StraightSkeleton2>(CGAL::create_interior_straight_skeleton_2(
      _poly.vertices_begin(), _poly.vertices_end()));
  }
};

typedef std::shared_ptr<Polygon2> Polygon2Ref;


class PolygonWithHoles2 {
  CGAL::Polygon_with_holes_2<K> _poly;

public:
  PolygonWithHoles2(const Polygon2Ref &p_poly) : _poly(p_poly->cgal_polygon_2()) {
  }

  void add_hole(const Polygon2Ref &p_poly) {
    _poly.add_hole(p_poly->cgal_polygon_2());
  }

  StraightSkeleton2Ref create_interior_straight_skeleton() const {
    return std::make_shared<StraightSkeleton2>(CGAL::create_interior_straight_skeleton_2(
      _poly.outer_boundary().vertices_begin(), _poly.outer_boundary().vertices_end(),
      _poly.holes_begin(), _poly.holes_end()));
  }
};

typedef std::shared_ptr<PolygonWithHoles2> PolygonWithHoles2Ref;



PYBIND11_MODULE(pycgalmini, m) {
  py::class_<Polygon2, Polygon2Ref> polygon2(m, "Polygon2D");
  polygon2.def(py::init<py::array_t<double>>());
  polygon2.def("create_interior_straight_skeleton", &Polygon2::create_interior_straight_skeleton);

  py::class_<PolygonWithHoles2, PolygonWithHoles2Ref> polygon_with_holes2(m, "PolygonWithHoles2D");
  polygon_with_holes2.def(py::init<Polygon2Ref>());
  polygon_with_holes2.def("create_interior_straight_skeleton", &PolygonWithHoles2::create_interior_straight_skeleton);

  py::class_<StraightSkeleton2, StraightSkeleton2Ref> skeleton2(m, "StraightSkeleton2D");
  skeleton2.def_property_readonly("bisectors", &StraightSkeleton2::bisectors);
}
