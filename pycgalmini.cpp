#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/create_straight_skeleton_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>

#include <CGAL/Boolean_set_operations_2.h>

namespace py = pybind11;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel KE;


class StraightSkeleton2 {
  const boost::shared_ptr<CGAL::Straight_skeleton_2<K>> _ref;

public:
  StraightSkeleton2(const boost::shared_ptr<CGAL::Straight_skeleton_2<K>> &ref) : _ref(ref) {
  }

  py::array_t<double> vertices() const {
    const size_t n = _ref->size_of_vertices();

    py::array_t<double> pts({ n, size_t(2) });
    auto r = pts.mutable_unchecked<2>();

    size_t j = 0;
    for (auto i = _ref->vertices_begin(); i != _ref->vertices_end(); i++, j++) {
      r(j, 0) = i->point().x();
      r(j, 1) = i->point().y();
    }

    return pts;
  }

  py::array_t<int32_t> bisectors() const {
    size_t n = 0;
    for (auto i = _ref->halfedges_begin(); i != _ref->halfedges_end(); i++) {
      if (i->is_bisector()) {
        n += 1;
      }
    }

    std::map<size_t, int32_t> vertex_ids;
    {
      size_t j = 0;
      for (auto i = _ref->vertices_begin(); i != _ref->vertices_end(); i++, j++) {
        vertex_ids[i->id()] = j;
      }    
    }

    py::array_t<int32_t> pts({ n, size_t(2) });
    auto r = pts.mutable_unchecked<2>();
    size_t j = 0;

    for (auto i = _ref->halfedges_begin(); i != _ref->halfedges_end(); i++) {
      if (i->is_bisector()) {
        r(j, 0) = vertex_ids[i->vertex()->id()];
        r(j, 1) = vertex_ids[i->opposite()->vertex()->id()];
        j += 1;
      }
    }

    assert(j == n);

    return pts;
  }
};

typedef std::shared_ptr<StraightSkeleton2> StraightSkeleton2Ref;


void _build_polygon(CGAL::Polygon_2<K> &p_poly, const py::array_t<double> &p_vertices) {
    auto r = p_vertices.unchecked<2>();
    if (r.shape(1) != 2) {
      throw std::runtime_error("vertices need to be 2 dimensional");
    }
    const ssize_t n = r.shape(0);

    p_poly.clear();
    for (ssize_t i = 0; i < n; i++) {
      p_poly.push_back(K::Point_2(r(i, 0), r(i, 1)));
    }
}


class Polygon2 {
  CGAL::Polygon_2<K> _poly;

public:
  inline Polygon2(const CGAL::Polygon_2<K> &p_poly) : _poly(p_poly) {
  }

  Polygon2(py::array_t<double> p_vertices) {
    _build_polygon(_poly, p_vertices);
  }

  inline const CGAL::Polygon_2<K> &cgal_polygon() const {
    return _poly;
  }

  py::array_t<double> coords() const {
    py::array_t<double> pts({ _poly.size(), size_t(2) });
    auto r = pts.mutable_unchecked<2>();
    size_t j = 0;

    for (auto i = _poly.vertices_begin(); i != _poly.vertices_end(); i++, j++) {
      r(j, 0) = i->x();
      r(j, 1) = i->y();
    }

    return pts;
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
  inline PolygonWithHoles2() {
  }

  inline PolygonWithHoles2(const CGAL::Polygon_with_holes_2<K> &p_poly) : _poly(p_poly) {
  }

  PolygonWithHoles2(py::array_t<double> p_vertices) {
    CGAL::Polygon_2<K> poly;
    _build_polygon(poly, p_vertices);
    _poly = CGAL::Polygon_with_holes_2<K>(poly);
  }

  inline CGAL::Polygon_with_holes_2<K> &mutable_cgal_polygon() {
    return _poly;
  }

  inline const CGAL::Polygon_with_holes_2<K> &cgal_polygon() const {
    return _poly;
  }

  Polygon2Ref exterior() const {
    return std::make_shared<Polygon2>(_poly.outer_boundary());
  }

  py::list interiors() const {
    py::list res;
    for (auto i = _poly.holes_begin(); i != _poly.holes_end(); i++) {
      res.append(std::make_shared<Polygon2>(*i));
    }
    return res;
  }

  void add_hole(py::array_t<double> p_vertices) {
    CGAL::Polygon_2<K> poly;
    _build_polygon(poly, p_vertices);
    _poly.add_hole(poly);
  }

  StraightSkeleton2Ref create_interior_straight_skeleton() const {
    return std::make_shared<StraightSkeleton2>(CGAL::create_interior_straight_skeleton_2(
      _poly.outer_boundary().vertices_begin(), _poly.outer_boundary().vertices_end(),
      _poly.holes_begin(), _poly.holes_end()));
  }
};

typedef std::shared_ptr<PolygonWithHoles2> PolygonWithHoles2Ref;


PolygonWithHoles2Ref random_polygon(int n, double extent) {
  CGAL::Polygon_2<K> poly;
  CGAL::random_polygon_2(n, std::back_inserter(poly), CGAL::Random_points_in_square_2<K::Point_2>(extent));

  py::array_t<double> pts({ poly.size(), size_t(2) });
  auto r = pts.mutable_unchecked<2>();
  size_t j = 0;
  for (auto i = poly.vertices_begin(); i != poly.vertices_end(); i++, j++) {
    r(j, 0) = i->x();
    r(j, 1) = i->y();
  }
  assert(j == poly.size());

  return std::make_shared<PolygonWithHoles2>(pts);
}

CGAL::Polygon_2<KE> to_exact(const CGAL::Polygon_2<K> &p_poly) {
  CGAL::Polygon_2<KE> exact;
  for (auto i = p_poly.vertices_begin(); i != p_poly.vertices_end(); i++) {
    exact.push_back(KE::Point_2(i->x(), i->y()));
  }
  return exact;
}

CGAL::Polygon_2<K> to_inexact(const CGAL::Polygon_2<KE> &p_poly) {
  CGAL::Polygon_2<K> inexact;
  for (auto i = p_poly.vertices_begin(); i != p_poly.vertices_end(); i++) {
    inexact.push_back(K::Point_2(i->x().exact().convert_to<double>(), i->y().exact().convert_to<double>()));
  }
  return inexact;
}

PolygonWithHoles2Ref join(PolygonWithHoles2Ref p_p1, PolygonWithHoles2Ref p_p2) {
  // https://doc.cgal.org/latest/Boolean_set_operations_2/index.html

  CGAL::Polygon_with_holes_2<KE> p1(to_exact(p_p1->cgal_polygon().outer_boundary()));
  CGAL::Polygon_with_holes_2<KE> p2(to_exact(p_p2->cgal_polygon().outer_boundary()));
  CGAL::Polygon_with_holes_2<KE> res;
  CGAL::join(p1, p2, res);

  CGAL::Polygon_with_holes_2<K> r(to_inexact(res.outer_boundary()));
  return std::make_shared<PolygonWithHoles2>(r);
}


PYBIND11_MODULE(pycgalmini, m) {
  py::class_<Polygon2, Polygon2Ref> polygon2(m, "LineString");
  polygon2.def(py::init<py::array_t<double>>());
  polygon2.def("create_interior_straight_skeleton", &Polygon2::create_interior_straight_skeleton);
  polygon2.def_property_readonly("coords", &Polygon2::coords);

  py::class_<PolygonWithHoles2, PolygonWithHoles2Ref> polygon_with_holes2(m, "Polygon");
  polygon_with_holes2.def(py::init<py::array_t<double>>());
  polygon_with_holes2.def_property_readonly("exterior", &PolygonWithHoles2::exterior);
  polygon_with_holes2.def_property_readonly("interiors", &PolygonWithHoles2::interiors);
  polygon_with_holes2.def("add_hole", &PolygonWithHoles2::add_hole);
  polygon_with_holes2.def("create_interior_straight_skeleton", &PolygonWithHoles2::create_interior_straight_skeleton);

  py::class_<StraightSkeleton2, StraightSkeleton2Ref> skeleton2(m, "StraightSkeleton");
  skeleton2.def_property_readonly("vertices", &StraightSkeleton2::vertices);
  skeleton2.def_property_readonly("bisectors", &StraightSkeleton2::bisectors);

  m.def("random_polygon", &random_polygon);

  m.def("join", &join);
}
