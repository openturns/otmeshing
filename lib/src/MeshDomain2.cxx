//                                               -*- C++ -*-
/**
 *  @brief Mesh domain with distance
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "otmeshing/MeshDomain2.hxx"

#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/SpecFunc.hxx>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#if CGAL_VERSION_NR >= 1060000000
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#else
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#endif
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(MeshDomain2)
static const Factory<MeshDomain2> Factory_MeshDomain2;


/* Default constructor */
MeshDomain2::MeshDomain2()
  : MeshDomain()
{
  // Nothing to do
}

/* Parameters constructor */
MeshDomain2::MeshDomain2(const OT::Mesh & mesh)
: MeshDomain(mesh)
{
  // Nothing to do
}

/* Virtual constructor */
MeshDomain2 * MeshDomain2::clone() const
{
  return new MeshDomain2(*this);
}

/* Compute the Euclidean distance from a given point to the domain */
Scalar MeshDomain2::computeDistance(const Point & point) const
{
  return computeDistance(Sample(1, point))(0, 0);
}

Sample MeshDomain2::computeDistance(const Sample & points) const
{
  const UnsignedInteger dimension = getDimension();
  const UnsignedInteger size = points.getSize();
  if (points.getDimension() != dimension)
    throw InvalidArgumentException(HERE) << "Expected a point of dimension " << dimension << " got " << points.getDimension();
  const Mesh mesh(getMesh());
  const Sample vertices(mesh.getVertices());
  const IndicesCollection simplices(mesh.getSimplices());
  Sample distances(size, 1);

  if (dimension == 2)
  {
    using KernelInexact = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point2 = KernelInexact::Point_2;
    using Point3 = KernelInexact::Point_3;
    using Segment2 = KernelInexact::Segment_2;
    using Segment3 = KernelInexact::Segment_3;
    using Ray2 = KernelInexact::Ray_2;
    using TriangleIterator = std::vector<Segment3>::iterator;
#if CGAL_VERSION_NR >= 1060000000
    using Primitive = CGAL::AABB_triangle_primitive_3<KernelInexact, TriangleIterator>;
    using AABB_Traits = CGAL::AABB_traits_3<KernelInexact, Primitive>;
#else
    using Primitive = CGAL::AABB_triangle_primitive<KernelInexact, TriangleIterator>;
    using AABB_Traits = CGAL::AABB_traits<KernelInexact, Primitive>;
#endif
    using Tree = CGAL::AABB_tree<AABB_Traits>;

    // we want to filter out internal facet
    // external facets are only referenced by one cell
    std::map<Indices, UnsignedInteger> facetMap;
    for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    {
      const UnsignedInteger i0 = simplices(i, 0);
      const UnsignedInteger i1 = simplices(i, 1);
      const UnsignedInteger i2 = simplices(i, 2);
      Indices f1 = {i0, i1};
      Indices f2 = {i1, i2};
      Indices f3 = {i2, i0};
      std::sort(f1.begin(), f1.end());
      std::sort(f2.begin(), f2.end());
      std::sort(f3.begin(), f3.end());
      ++ facetMap[f1];
      ++ facetMap[f2];
      ++ facetMap[f3];
    }

    // Build list of boundary edges
    std::vector<Segment2> boundary_edges2;
    std::vector<Segment3> boundary_edges3;
    for (const auto & elt : facetMap)
    {
      if (elt.second == 1) // boundary edge
      {
        const UnsignedInteger i0 = elt.first[0];
        const UnsignedInteger i1 = elt.first[1];

        const Point2 v0{vertices(i0, 0), vertices(i0, 1)};
        const Point2 v1{vertices(i1, 0), vertices(i1, 1)};

        boundary_edges2.emplace_back(v0, v1);
        boundary_edges3.emplace_back(Point3(v0[0], v0[1], 0.0), Point3(v1[0], v1[1], 0.0));
      }
    }

    // Build 3d tree (CGAL<6 does not support 2d AABB tree, so lift points)
    Tree tree(boundary_edges3.begin(), boundary_edges3.end());
    tree.accelerate_distance_queries();

    for (UnsignedInteger i = 0; i < size; ++ i)
    {
      // distance to closest simplex
      const Point3 query(points(i, 0), points(i, 1), 0.0);
      const Point3 closest = tree.closest_point(query);
      distances(i, 0) = std::sqrt(CGAL::squared_distance(query, closest));

      // Ray casting for general boundary edge list
      UnsignedInteger intersections = 0;
      const Ray2 ray(Point2(points(i, 0), points(i, 1)),
                     Point2(points(i, 0) + 1.0, points(i, 1))); // horizontal ray
      for (const auto& edge2 : boundary_edges2)
      {
        if (CGAL::do_intersect(ray, edge2))
          ++ intersections;
      }
      const Bool inside = (intersections % 2) == 1;
      if (inside)
        distances(i, 0) = -distances(i, 0);
      LOGDEBUG(OSS() << "query=" << points[i] << " closest=" << Point({closest[0], closest[1]}) << " intersections=" << intersections << " inside=" << inside << " dist=" << distances(i, 0));
    }
  }
  else if (dimension == 3)
  {
    using KernelInexact = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point3 = KernelInexact::Point_3;
    using Mesh3 = CGAL::Surface_mesh<Point3>;
    using Side_of_triangle_mesh = CGAL::Side_of_triangle_mesh<Mesh3, KernelInexact>;
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh3>;
#if CGAL_VERSION_NR >= 1060000000
    using AABB_Traits = CGAL::AABB_traits_3<KernelInexact, Primitive>;
#else
    using AABB_Traits = CGAL::AABB_traits<KernelInexact, Primitive>;
#endif
    using Tree = CGAL::AABB_tree<AABB_Traits>;

    // we want to filter out internal facet
    // external facets are only referenced by one cell
    std::map<Indices, UnsignedInteger> facetMap;
    for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    {
      const UnsignedInteger i0 = simplices(i, 0);
      const UnsignedInteger i1 = simplices(i, 1);
      const UnsignedInteger i2 = simplices(i, 2);
      const UnsignedInteger i3 = simplices(i, 3);
      Indices f1 = {i0, i1, i2};
      Indices f2 = {i0, i1, i3};
      Indices f3 = {i0, i2, i3};
      Indices f4 = {i1, i2, i3};
      std::sort(f1.begin(), f1.end());
      std::sort(f2.begin(), f2.end());
      std::sort(f3.begin(), f3.end());
      std::sort(f4.begin(), f4.end());
      ++ facetMap[f1];
      ++ facetMap[f2];
      ++ facetMap[f3];
      ++ facetMap[f4];
    }

    // build boundary mesh representation instead of triangle representation
    // allows to use Side_of_triangle_mesh boundary point location
    Mesh3 mesh3;
    std::map<Point3, Mesh3::Vertex_index> vMap;
    for (const auto & elt : facetMap)
    {
      // check if external facet
      if (elt.second == 1)
      {
        const UnsignedInteger i0 = elt.first[0];
        const UnsignedInteger i1 = elt.first[1];
        const UnsignedInteger i2 = elt.first[2];

        const Point3 v0{vertices(i0, 0), vertices(i0, 1), vertices(i0, 2)};
        const Point3 v1{vertices(i1, 0), vertices(i1, 1), vertices(i1, 2)};
        const Point3 v2{vertices(i2, 0), vertices(i2, 1), vertices(i2, 2)};

        // only add vertices of boundary facets
        if (vMap.find(v0) == vMap.end())
          vMap[v0] = mesh3.add_vertex(v0);
        if (vMap.find(v1) == vMap.end())
          vMap[v1] = mesh3.add_vertex(v1);
        if (vMap.find(v2) == vMap.end())
          vMap[v2] = mesh3.add_vertex(v2);

        Mesh3::Vertex_index v[3] = {vMap[v0], vMap[v1], vMap[v2]};
        auto f = mesh3.add_face(v[0], v[1], v[2]);

        // check for non-manifold edge
        if (f == Mesh3::null_face())
          mesh3.add_face(v[0], v[2], v[1]);
      }
    }
    if (!CGAL::is_closed(mesh3))
      throw InternalException(HERE) << "MeshDomain2.computeDistance(3d): mesh should be closed";

    // build tree
    Tree tree(mesh3.faces().begin(), mesh3.faces().end(), mesh3);
    tree.accelerate_distance_queries();

    // this is more robust than using simple ray intersection with AABB_tree
    const Side_of_triangle_mesh insideTester(mesh3);

    for (UnsignedInteger i = 0; i < size; ++ i)
    {
      // distance to closest facet
      const Point3 query(points(i, 0), points(i, 1), points(i, 2));
      const Point3 closest = tree.closest_point(query);
      distances(i, 0) = std::sqrt(CGAL::squared_distance(query, closest));

      // check whether the point is inside/outside the mesh
      const CGAL::Bounded_side side = insideTester(query);
      distances(i, 0) = (side == CGAL::ON_BOUNDED_SIDE) ? -distances(i, 0) : distances(i, 0);

      LOGDEBUG(OSS() << "query=" << points[i] << " closest=" << Point({closest[0], closest[1], closest[2]}) << " inside=" << (side == CGAL::ON_BOUNDED_SIDE));
    }
  }
  else
    throw NotYetImplementedException(HERE) << "MeshDomain2 does not support dimension " << dimension;
  return distances;
}

}
