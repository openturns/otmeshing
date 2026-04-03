//                                               -*- C++ -*-
/**
 *  @brief Meshing algorithm for points
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "otmeshing/ConvexDecompositionMesher.hxx"
#include "otmeshing/CloudMesher.hxx"

#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/SpecFunc.hxx>

#include <iterator>
#include <map>
#include <unordered_map>

// 3d
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// 2d
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>

#include <queue>
#include <set>

#ifdef OPENTURNS_HAVE_COACD
#include <CoACD/coacd.h>
#endif

using namespace OT;



namespace OTMESHING
{

CLASSNAMEINIT(ConvexDecompositionMesher)

static Factory<ConvexDecompositionMesher> Factory_ConvexDecompositionMesher;


/* Default constructor */
ConvexDecompositionMesher::ConvexDecompositionMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor method */
ConvexDecompositionMesher * ConvexDecompositionMesher::clone() const
{
  return new ConvexDecompositionMesher(*this);
}


Collection<Mesh> ConvexDecompositionMesher::build(const Mesh & mesh) const
{
  const UnsignedInteger dimension = mesh.getDimension();
  const UnsignedInteger intrinsicDimension = mesh.getIntrinsicDimension();
  const Sample vertices(mesh.getVertices());
  const IndicesCollection simplices(mesh.getSimplices());
  Collection<Mesh> result;

  if (dimension == 2 && !useSimplicesDecomposition_)
  {
    using KernelInexact = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point_2 = KernelInexact::Point_2;
    using Polygon_2 = CGAL::Polygon_2<KernelInexact>;
    using Point3 = KernelInexact::Point_3;
    using Mesh3 = CGAL::Surface_mesh<Point3>;

    Mesh3 mesh3;
    std::map<Point3, Mesh3::Vertex_index> vMap;
    for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    {
      const UnsignedInteger i0 = simplices(i, 0);
      const UnsignedInteger i1 = simplices(i, 1);
      const UnsignedInteger i2 = simplices(i, 2);
      const Point3 v0{vertices(i0, 0), vertices(i0, 1), 0.0};
      const Point3 v1{vertices(i1, 0), vertices(i1, 1), 0.0};
      const Point3 v2{vertices(i2, 0), vertices(i2, 1), 0.0};

      if (vMap.find(v0) == vMap.end())
        vMap[v0] = mesh3.add_vertex(v0);
      if (vMap.find(v1) == vMap.end())
        vMap[v1] = mesh3.add_vertex(v1);
      if (vMap.find(v2) == vMap.end())
        vMap[v2] = mesh3.add_vertex(v2);

      if (mesh3.add_face(vMap[v0], vMap[v1], vMap[v2]) == Mesh3::null_face())
        throw InternalException(HERE) << "Degenerate or duplicate face";
    }

    // decompose into connected components
    std::vector<std::size_t> components(mesh3.number_of_faces());
    const UnsignedInteger componentsNumber = CGAL::Polygon_mesh_processing::connected_components(mesh3, CGAL::make_property_map(components));
    LOGDEBUG(OSS() << "Number of connected components" << componentsNumber);

    std::vector<Point_2> vertices2;
    for (UnsignedInteger i = 0; i < vertices.getSize(); ++ i)
    {
      const Point_2 v0{vertices(i, 0), vertices(i, 1)};
      vertices2.push_back(v0);
    }

    // for each disconnected components
    CloudMesher mesher;
    for (UnsignedInteger cid = 0; cid < componentsNumber; ++ cid)
    {
      // collect triangle indices for this component
      Indices componentTriangles;
      for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
        if (components[Mesh3::Face_index(i)] == cid)
          componentTriangles.add(i);

      // greedy algorithm to form convex components starting from each triangle
      Indices used(simplices.getSize());
      for (const UnsignedInteger idx : componentTriangles)
      {
        if (used[idx]) continue;

        // initial triangle to start the convex polygon from
        Polygon_2 poly;
        poly.push_back(vertices2[simplices(idx, 0)]);
        poly.push_back(vertices2[simplices(idx, 1)]);
        poly.push_back(vertices2[simplices(idx, 2)]);
        used[idx] = 1;

        // continue while we succedeed in aggregating a triangle
        Bool mergedAny = false;
        do
        {
          mergedAny = false;
          for (const UnsignedInteger jdx : componentTriangles)
          {
            if (used[jdx]) continue;

            // check if triangle shares an edge with the current polygon
            // by checking if two of its vertices correspond to consecutive polygon vertices
            const Point_2 q0{vertices2[simplices(jdx, 0)]};
            const Point_2 q1{vertices2[simplices(jdx, 1)]};
            const Point_2 q2{vertices2[simplices(jdx, 2)]};

            for (auto vi = poly.vertices_begin(); vi != poly.vertices_end(); ++ vi)
            {
              auto next = std::next(vi);
              if (next == poly.vertices_end())
                next = poly.vertices_begin();

              Point_2 newPoint;
              Bool match = false;
              if ((*vi == q0 && *next == q1) || (*vi == q1 && *next == q0)) { newPoint = q2; match = true; }
              else if ((*vi == q1 && *next == q2) || (*vi == q2 && *next == q1)) { newPoint = q0; match = true; }
              else if ((*vi == q2 && *next == q0) || (*vi == q0 && *next == q2)) { newPoint = q1; match = true; }

              if (match)
              {
                Polygon_2 merged{poly};
                // find the same edge in the copy and insert the new vertex
                for (auto mvi = merged.vertices_begin(); mvi != merged.vertices_end(); ++ mvi)
                {
                  auto mnext = std::next(mvi);
                  if (mnext == merged.vertices_end())
                    mnext = merged.vertices_begin();
                  if ((*mvi == *vi && *mnext == *next) || (*mvi == *next && *mnext == *vi))
                  {
                    merged.insert(mnext, newPoint);
                    break;
                  }
                }
                if (merged.is_convex())
                {
                  poly = merged;
                  used[jdx] = 1;
                  mergedAny = true;
                }
                break;
              }
            }
          }
        } while (mergedAny);

        Sample verticesI(0, dimension);
        for (const auto & p : poly)
          verticesI.add(Point({p[0], p[1]}));
        result.add(mesher.build(verticesI));
      }
    }
  }
  else if (dimension == 3 && !useSimplicesDecomposition_)
  {
    if (intrinsicDimension == 2)
    {
      // Surface mesh: use CGAL exact convex decomposition (coacd returns open patches)
      using KernelExact = CGAL::Exact_predicates_exact_constructions_kernel;
      using Polyhedron3 = CGAL::Polyhedron_3<KernelExact>;
      using HDS = Polyhedron3::HalfedgeDS;
      using Nef_polyhedron3 = CGAL::Nef_polyhedron_3<KernelExact>;
      using Point_3 = KernelExact::Point_3;

      // Build Nef polyhedron from surface mesh
      Polyhedron3 poly;
      CGAL::Polyhedron_incremental_builder_3<HDS> builder(poly.hds(), true);
      builder.begin_surface(vertices.getSize(), simplices.getSize());
      for (UnsignedInteger i = 0; i < vertices.getSize(); ++ i)
      {
        const Point_3 p{vertices(i, 0), vertices(i, 1), vertices(i, 2)};
        builder.add_vertex(p);
      }
      for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
      {
        builder.begin_facet();
        for (UnsignedInteger j = 0; j < dimension; ++ j)
          builder.add_vertex_to_facet(simplices(i, j));
        builder.end_facet();
      }
      builder.end_surface();
      if (builder.check_unconnected_vertices())
      {
        LOGINFO("ConvexDecompositionMesher detected unconnected vertices, removing");
        if (!builder.remove_unconnected_vertices())
          throw InvalidArgumentException(HERE) << "Polyhedron could not remove all unconnected vertices";
      }
      Nef_polyhedron3 nef(poly);

      // Extract convex components
      CGAL::convex_decomposition_3(nef);

      // the first volume is the outer volume, which is ignored in the decomposition
      for (auto ci = ++nef.volumes_begin(); ci != nef.volumes_end(); ++ci)
      {
        if (ci->mark())
        {
          Polyhedron3 part;
          nef.convert_inner_shell_to_polyhedron(ci->shells_begin(), part);
          if (part.empty()) continue;

          CGAL::Polygon_mesh_processing::triangulate_faces(part);

          Sample verticesI(part.size_of_vertices(), dimension);
          std::unordered_map<typename Polyhedron3::Vertex_const_handle, UnsignedInteger> vertexToIndexMap;
          UnsignedInteger vertexIndex = 0;
          for (auto vi = part.vertices_begin(); vi != part.vertices_end(); ++ vi)
          {
            const Point_3 & p = vi->point();
            for (UnsignedInteger j = 0; j < dimension; ++ j)
              verticesI(vertexIndex, j) = CGAL::to_double(p[j]);
            vertexToIndexMap[vi] = vertexIndex;
            ++ vertexIndex;
          }

          // arbitrarily select the apex as the first vertex
          const UnsignedInteger apexIndex = vertexToIndexMap[part.vertices_begin()];

          // build simplices
          Collection<Indices> simplexColl;
          for (auto f = part.facets_begin(); f != part.facets_end(); ++f)
          {
            auto h = f->facet_begin();
            Bool ok = true;
            for (UnsignedInteger j = 0; j < dimension; ++ j)
            {
              if (vertexToIndexMap[h->vertex()] == apexIndex)
              {
                ok = false;
                break;
              }
              ++ h;
            }
            if (ok)
            {
              h = f->facet_begin();
              Indices simplex(dimension + 1);
              simplex[0] = apexIndex;
              for (UnsignedInteger j = 0; j < dimension; ++ j)
              {
                simplex[j + 1] = vertexToIndexMap[h->vertex()];
                ++ h;
              }
              simplexColl.add(simplex);
            }
          }

          // remove unused vertices
          Indices usedVertices(verticesI.getSize());
          for (const auto & simplex : simplexColl)
            for (const UnsignedInteger vi : simplex)
              usedVertices[vi] = 1;

          Indices oldToNew(verticesI.getSize());
          Sample verticesICompact(0, dimension);
          for (UnsignedInteger i = 0; i < verticesI.getSize(); ++ i)
          {
            if (usedVertices[i])
            {
              Point p(dimension);
              for (UnsignedInteger j = 0; j < dimension; ++ j)
                p[j] = verticesI(i, j);
              verticesICompact.add(p);
              oldToNew[i] = verticesICompact.getSize() - 1;
            }
          }

          Collection<Indices> simplexCollCompact;
          for (const auto & simplex : simplexColl)
          {
            Indices simplexNew(simplex.getSize());
            for (UnsignedInteger j = 0; j < simplex.getSize(); ++ j)
              simplexNew[j] = oldToNew[simplex[j]];
            simplexCollCompact.add(simplexNew);
          }
          result.add(Mesh(verticesICompact, IndicesCollection(simplexCollCompact)));
        }
      } // for nef.volumes
      return result;
    }
    else if (intrinsicDimension == 3)
    {
#ifdef OPENTURNS_HAVE_COACD
      // tetra mesh -> extract external facets
      // external facets are only referenced by one cell
      // First pass: count using sorted keys (orientation-independent)
      const Point simplicesVolume(mesh.computeSimplicesVolume());
      std::map<Indices, UnsignedInteger> facetMap;
      for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
      {
        if (simplicesVolume[i] <= 0.0)
          continue;
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

      // Second pass: collect external triangles with outward-facing orientation
      struct Tri { int v[3]; };
      std::vector<Tri> triangles;
      for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
      {
        if (simplicesVolume[i] <= 0.0)
          continue;
        const UnsignedInteger i0 = simplices(i, 0);
        const UnsignedInteger i1 = simplices(i, 1);
        const UnsignedInteger i2 = simplices(i, 2);
        const UnsignedInteger i3 = simplices(i, 3);
        // Each face paired with its opposite vertex
        Indices faceIndices[4] = {{i0, i1, i2}, {i0, i1, i3}, {i0, i2, i3}, {i1, i2, i3}};
        UnsignedInteger opp[4] = {i3, i2, i1, i0};
        for (UnsignedInteger f = 0; f < 4; ++ f)
        {
          auto & face = faceIndices[f];
          Indices key = {face[0], face[1], face[2]};
          std::sort(key.begin(), key.end());
          if (facetMap[key] != 1) continue;

          // Orient face so normal points away from the tetrahedron interior
          // Compute face normal from current vertex order
          Scalar x0 = vertices(face[0], 0), y0 = vertices(face[0], 1), z0 = vertices(face[0], 2);
          Scalar x1 = vertices(face[1], 0), y1 = vertices(face[1], 1), z1 = vertices(face[1], 2);
          Scalar x2 = vertices(face[2], 0), y2 = vertices(face[2], 1), z2 = vertices(face[2], 2);
          Scalar nx = (y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0);
          Scalar ny = (z1 - z0) * (x2 - x0) - (x1 - x0) * (z2 - z0);
          Scalar nz = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
          // Vector from face center to opposite vertex
          Scalar cx = (x0 + x1 + x2) / 3.0;
          Scalar cy = (y0 + y1 + y2) / 3.0;
          Scalar cz = (z0 + z1 + z2) / 3.0;
          Scalar dx = vertices(opp[f], 0) - cx;
          Scalar dy = vertices(opp[f], 1) - cy;
          Scalar dz = vertices(opp[f], 2) - cz;
          // If normal points toward opposite vertex, flip it
          if (nx * dx + ny * dy + nz * dz > 0.0)
            std::swap(face[1], face[2]);

          Tri t;
          t.v[0] = static_cast<int>(face[0]);
          t.v[1] = static_cast<int>(face[1]);
          t.v[2] = static_cast<int>(face[2]);
          triangles.push_back(t);
        }
      }

      if (triangles.empty())
        throw InternalException(HERE) << "ConvexDecompositionMesher: no external facets found";

      // group triangles into connected components (by shared edges)
      std::map<std::pair<int, int>, std::vector<int>> edgeMap;
      for (int i = 0; i < static_cast<int>(triangles.size()); ++ i)
      {
        const auto & t = triangles[i];
        std::pair<int, int> edges[3] = {
          {std::min(t.v[0], t.v[1]), std::max(t.v[0], t.v[1])},
          {std::min(t.v[1], t.v[2]), std::max(t.v[1], t.v[2])},
          {std::min(t.v[0], t.v[2]), std::max(t.v[0], t.v[2])}
        };
        for (auto & e : edges)
          edgeMap[e].push_back(i);
      }

      std::vector<int> component(triangles.size(), -1);
      int componentCount = 0;
      for (int i = 0; i < static_cast<int>(triangles.size()); ++ i)
      {
        if (component[i] >= 0) continue;
        std::queue<int> q;
        q.push(i);
        component[i] = componentCount;
        while (!q.empty())
        {
          const int t = q.front(); q.pop();
          const auto & tri = triangles[t];
          std::pair<int, int> edges[3] = {
            {std::min(tri.v[0], tri.v[1]), std::max(tri.v[0], tri.v[1])},
            {std::min(tri.v[1], tri.v[2]), std::max(tri.v[1], tri.v[2])},
            {std::min(tri.v[0], tri.v[2]), std::max(tri.v[0], tri.v[2])}
          };
          for (auto & e : edges)
          {
            for (int adj : edgeMap[e])
              if (component[adj] < 0)
              {
                component[adj] = componentCount;
                q.push(adj);
              }
          }
        }
        ++ componentCount;
      }
      LOGDEBUG(OSS() << "External facet connected components: " << componentCount);

      // run coacd on each connected component separately
      for (int cid = 0; cid < componentCount; ++ cid)
      {
        // collect vertices referenced by this component
        std::set<int> vertSet;
        std::vector<Tri> compTriangles;
        for (int i = 0; i < static_cast<int>(triangles.size()); ++ i)
        {
          if (component[i] != cid) continue;
          const auto & t = triangles[i];
          vertSet.insert(t.v[0]);
          vertSet.insert(t.v[1]);
          vertSet.insert(t.v[2]);
          compTriangles.push_back(t);
        }

        if (compTriangles.empty()) continue;

        // build a local vertex index map
        std::map<int, int> globalToLocal;
        std::vector<int> localToGlobal;
        for (int gv : vertSet)
        {
          globalToLocal[gv] = localToGlobal.size();
          localToGlobal.push_back(gv);
        }

        coacd::Mesh input;
        input.vertices.resize(localToGlobal.size());
        for (std::size_t i = 0; i < localToGlobal.size(); ++ i)
        {
          input.vertices[i][0] = vertices(localToGlobal[i], 0);
          input.vertices[i][1] = vertices(localToGlobal[i], 1);
          input.vertices[i][2] = vertices(localToGlobal[i], 2);
        }
        for (auto & t : compTriangles)
        {
          input.indices.push_back({globalToLocal[t.v[0]], globalToLocal[t.v[1]], globalToLocal[t.v[2]]});
        }

        coacd::set_log_level(Log::HasDebug() ? "debug" : "off");
        const Scalar threshold = ResourceMap::GetAsScalar("ConvexDecompositionMesher-Threshold");
        std::vector<coacd::Mesh> output = coacd::CoACD(input, threshold);
        LOGDEBUG(OSS() << "Component " << cid << " N CONVEX=" << output.size());

        for (const auto & part : output)
        {
          // tetrahedralize convex piece via fan from centroid
          const UnsignedInteger nv = part.vertices.size();
          Scalar centroidX = 0.0, centroidY = 0.0, centroidZ = 0.0;
          Sample verticesI(nv + 1, dimension);
          for (UnsignedInteger i = 0; i < nv; ++ i)
          {
            verticesI(i, 0) = part.vertices[i][0];
            verticesI(i, 1) = part.vertices[i][1];
            verticesI(i, 2) = part.vertices[i][2];
            centroidX += part.vertices[i][0];
            centroidY += part.vertices[i][1];
            centroidZ += part.vertices[i][2];
          }
          centroidX /= nv;
          centroidY /= nv;
          centroidZ /= nv;
          verticesI(nv, 0) = centroidX;
          verticesI(nv, 1) = centroidY;
          verticesI(nv, 2) = centroidZ;
          const UnsignedInteger apexIndex = nv;
          Collection<Indices> simplexColl;
          for (const auto & tri : part.indices)
          {
            Indices simplex(4);
            simplex[0] = apexIndex;
            simplex[1] = static_cast<UnsignedInteger>(tri[0]);
            simplex[2] = static_cast<UnsignedInteger>(tri[1]);
            simplex[3] = static_cast<UnsignedInteger>(tri[2]);
            // Ensure positive signed volume
            const Scalar x0 = verticesI(simplex[1], 0) - verticesI(simplex[0], 0);
            const Scalar y0 = verticesI(simplex[1], 1) - verticesI(simplex[0], 1);
            const Scalar z0 = verticesI(simplex[1], 2) - verticesI(simplex[0], 2);
            const Scalar x1 = verticesI(simplex[2], 0) - verticesI(simplex[0], 0);
            const Scalar y1 = verticesI(simplex[2], 1) - verticesI(simplex[0], 1);
            const Scalar z1 = verticesI(simplex[2], 2) - verticesI(simplex[0], 2);
            const Scalar x2 = verticesI(simplex[3], 0) - verticesI(simplex[0], 0);
            const Scalar y2 = verticesI(simplex[3], 1) - verticesI(simplex[0], 1);
            const Scalar z2 = verticesI(simplex[3], 2) - verticesI(simplex[0], 2);
            if (x0 * (y1 * z2 - z1 * y2) - y0 * (x1 * z2 - z1 * x2) + z0 * (x1 * y2 - y1 * x2) < 0.0)
              std::swap(simplex[2], simplex[3]);
            simplexColl.add(simplex);
          }
          if (!simplexColl.isEmpty())
            result.add(Mesh(verticesI, IndicesCollection(simplexColl)));
        }
      } // for components

#else // !OPENTURNS_HAVE_COACD -> CGAL exact decomposition for volumetric meshes
    // build from the volumetric mesh
    using KernelExact = CGAL::Exact_predicates_exact_constructions_kernel;
    using Polyhedron = CGAL::Polyhedron_3<KernelExact>;
    using Nef_polyhedron = CGAL::Nef_polyhedron_3<KernelExact>;
    using Point_3 = KernelExact::Point_3;

    Nef_polyhedron nef;
    const Point simplicesVolume(mesh.computeSimplicesVolume());
    for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    {
      if (simplicesVolume[i] <= 0.0)
        continue;

      const UnsignedInteger i0 = simplices(i, 0);
      const UnsignedInteger i1 = simplices(i, 1);
      const UnsignedInteger i2 = simplices(i, 2);
      const UnsignedInteger i3 = simplices(i, 3);

      const Point_3 v0{vertices(i0, 0), vertices(i0, 1), vertices(i0, 2)};
      const Point_3 v1{vertices(i1, 0), vertices(i1, 1), vertices(i1, 2)};
      const Point_3 v2{vertices(i2, 0), vertices(i2, 1), vertices(i2, 2)};
      const Point_3 v3{vertices(i3, 0), vertices(i3, 1), vertices(i3, 2)};

      Polyhedron poly;
      poly.make_tetrahedron(v0, v1, v2, v3);
      const Nef_polyhedron tetra(poly);
      nef += tetra;
    }

    // Extract convex components
    CGAL::convex_decomposition_3(nef);

    // the first volume is the outer volume, which is ignored in the decomposition
    for (auto ci = ++nef.volumes_begin(); ci != nef.volumes_end(); ++ci)
    {
      if (ci->mark())
      {
        Polyhedron part;
        nef.convert_inner_shell_to_polyhedron(ci->shells_begin(), part);

        if (part.empty())
          continue;

        CGAL::Polygon_mesh_processing::triangulate_faces(part);

        Sample verticesI(part.size_of_vertices(), dimension);
        std::unordered_map<typename Polyhedron::Vertex_const_handle, UnsignedInteger> vertexToIndexMap;
        UnsignedInteger vertexIndex = 0;
        for (auto vi = part.vertices_begin(); vi != part.vertices_end(); ++ vi)
        {
          const Point_3 & p = vi->point();
          for (UnsignedInteger j = 0; j < dimension; ++ j)
            verticesI(vertexIndex, j) = CGAL::to_double(p[j]);
          vertexToIndexMap[vi] = vertexIndex;
          ++ vertexIndex;
        }

        // arbitrarily select the apex as the first vertex
        const UnsignedInteger apexIndex = vertexToIndexMap[part.vertices_begin()];

        // build simplices
        Collection<Indices> simplexColl;
        for (auto f = part.facets_begin(); f != part.facets_end(); ++f)
        {
          auto h = f->facet_begin();

          // filter out the facets incident to the apex vertex
          Bool ok = true;
          for (UnsignedInteger j = 0; j < dimension; ++ j)
          {
            if (vertexToIndexMap[h->vertex()] == apexIndex)
            {
              ok = false;
              break;
            }
            ++ h;
          }

          if (ok)
          {
            h = f->facet_begin();
            Indices simplex(dimension + 1);
            simplex[0] = apexIndex;
            for (UnsignedInteger j = 0; j < dimension; ++ j)
            {
              simplex[j + 1] = vertexToIndexMap[h->vertex()];
              ++ h;
            }
            simplexColl.add(simplex);
          }
        }

        // remove unused vertices (some may only appear in facets incident to the apex)
        Indices usedVertices(verticesI.getSize());
        for (const auto & simplex : simplexColl)
          for (const UnsignedInteger vi : simplex)
            usedVertices[vi] = 1;

        Indices oldToNew(verticesI.getSize());
        Sample verticesICompact(0, dimension);
        for (UnsignedInteger i = 0; i < verticesI.getSize(); ++ i)
        {
          if (usedVertices[i])
          {
            Point p(dimension);
            for (UnsignedInteger j = 0; j < dimension; ++ j)
              p[j] = verticesI(i, j);
            verticesICompact.add(p);
            oldToNew[i] = verticesICompact.getSize() - 1;
          }
        }

        Collection<Indices> simplexCollCompact;
        for (const auto & simplex : simplexColl)
        {
          Indices simplexNew(simplex.getSize());
          for (UnsignedInteger j = 0; j < simplex.getSize(); ++ j)
            simplexNew[j] = oldToNew[simplex[j]];
          simplexCollCompact.add(simplexNew);
        }
        result.add(Mesh(verticesICompact, IndicesCollection(simplexCollCompact)));
      } // if mark
    } // for nef
#endif // OPENTURNS_HAVE_COACD
    } // intrinsic dim=3
    else
      throw InvalidArgumentException(HERE) << "ConvexDecompositionMesher expected intrinsic dimension=2|3 got " << intrinsicDimension;
  } // dim = 3
  else if (dimension == intrinsicDimension)
  {
    const Point simplicesVolume(mesh.computeSimplicesVolume());

    // LevelSetMesher can yield almost empty cells
    // possible workaround with key LevelSetMesher-SolveEquation=False
    const Scalar smallVolume = simplicesVolume.norm1() * SpecFunc::Precision;

    Indices simplex(dimension + 1);
    Indices standardSimplex(dimension + 1);
    standardSimplex.fill();
    const IndicesCollection uniqueSimplex(1, dimension + 1, standardSimplex);
    Sample simplexVertices(dimension + 1, dimension);
    for (UnsignedInteger simplexIndex = 0; simplexIndex < simplices.getSize(); ++ simplexIndex)
    {
      // Skip small simplices
      if (!(simplicesVolume[simplexIndex] > smallVolume))
        continue;
      std::copy(simplices.cbegin_at(simplexIndex), simplices.cend_at(simplexIndex), simplex.begin());

      // Here we should keep the vertices untouched in order to avoid partial copies
      // but unused vertices throw an exception in the validity check
      for (UnsignedInteger j = 0; j <= dimension; ++j)
      {
        const UnsignedInteger localJ = simplex[j];
        for (UnsignedInteger k = 0; k < dimension; ++k)
          simplexVertices(j, k) = vertices(localJ, k);
      }
      result.add(Mesh(simplexVertices, uniqueSimplex));
    }
  }
  else
    throw InvalidArgumentException(HERE) << "ConvexDecompositionMesher expected dimension=3 and intrinsicDimension = 2|3, or dimension=intrinsicDimension, here got dimension=" << dimension << " and intrinsicDimension=" << intrinsicDimension;
  return result;
}

/* Check if mesh is convex */
Bool ConvexDecompositionMesher::IsConvex(const Mesh & mesh)
{
  CloudMesher mesher;
  const Scalar vm = mesh.getVolume();
  const Scalar vc = mesher.build(mesh.getVertices()).getVolume();
  return (vc > 0.0) && (std::abs((vm - vc) / vc) < std::sqrt(SpecFunc::Precision));
}

/* Simplices decomposition flag */
void ConvexDecompositionMesher::setUseSimplicesDecomposition(const Bool useSimplicesDecomposition)
{
  useSimplicesDecomposition_ = useSimplicesDecomposition;
}

Bool ConvexDecompositionMesher::getUseSimplicesDecomposition() const
{
  return useSimplicesDecomposition_;
}

/* String converter */
String ConvexDecompositionMesher::__repr__() const
{
  OSS oss;
  oss << "class=" << ConvexDecompositionMesher::GetClassName();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void ConvexDecompositionMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("useSimplicesDecomposition_", useSimplicesDecomposition_);
}

/* Method load() reloads the object from the StorageManager */
void ConvexDecompositionMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("useSimplicesDecomposition_", useSimplicesDecomposition_);
}


} /* namespace OTMESHING */
