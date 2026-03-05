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

  if (dimension == 2)
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

            // check if triangles share and edge
            Indices sharedVertexJIndex;
            for (UnsignedInteger i1 = 0; i1 < 3; ++ i1)
              for (UnsignedInteger i2 = 0; i2 < 3; ++ i2)
                if (simplices(idx, i1) == simplices(jdx, i2))
                  sharedVertexJIndex.add(i2);
            const Bool shareEdge = sharedVertexJIndex.getSize() == 2;
            if (!shareEdge) continue;

            // merge triangle jdx by inserting the vertex that is not shared in the shared edge
            const UnsignedInteger newVertexIndex = simplices(jdx, sharedVertexJIndex.complement(3)[0]);
            const Indices commonVertexIndex = {simplices(jdx, sharedVertexJIndex[0]), simplices(jdx, sharedVertexJIndex[1])};
            const Point_2 p0{vertices2[commonVertexIndex[0]]};
            const Point_2 p1{vertices2[commonVertexIndex[1]]};
            const Point_2 newPoint{vertices2[newVertexIndex]};
            Polygon_2 merged{poly};
            for (auto vi = merged.vertices_begin(); vi != merged.vertices_end(); ++ vi)
            {
              // wrap around
              auto next = std::next(vi);
              if (next == merged.vertices_end())
                next = merged.vertices_begin();

              if ((*vi == p0 && *next == p1) || (*vi == p1 && *next == p0))
              {
                merged.insert(next, newPoint);
                break;
              }
            }
            if (merged.is_convex())
            {
              poly = merged;
              used[jdx] = 1;
              mergedAny = true;
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
  else if (dimension == 3)
  {
    using KernelExact = CGAL::Exact_predicates_exact_constructions_kernel;
    using Polyhedron = CGAL::Polyhedron_3<KernelExact>;
    using HDS = Polyhedron::HalfedgeDS;
    using Nef_polyhedron = CGAL::Nef_polyhedron_3<KernelExact>;
    using Point_3 = KernelExact::Point_3;

    Nef_polyhedron nef;
    if (intrinsicDimension == 2)
    {
      // build from the surface mesh
      Polyhedron poly;
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

      // we have to check unconnected vertices otherwise older CGAL crashes when converting to Nef_polyhedron
      if (builder.check_unconnected_vertices())
      {
        LOGINFO("ConvexDecompositionMesher detected unconnected vertices, removing");
        if (!builder.remove_unconnected_vertices())
          throw InvalidArgumentException(HERE) << "Polyhedron could not remove all unconnected vertices";
      }

      if (!poly.is_valid(Log::HasDebug()))
        throw InvalidArgumentException(HERE) << "Polyhedron must be valid";

      if (!poly.is_closed())
        throw InvalidArgumentException(HERE) << "Polyhedron must be closed";

      nef = Nef_polyhedron(poly);
    }
    else if (intrinsicDimension == 3)
    {
      // build from the volumetric mesh
      for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
      {
        // cannot exclude small valume tetras as the nef can loose its 2-manifold property

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
    }
    else
      throw InvalidArgumentException(HERE) << "ConvexDecompositionMesher expected intrinsic dimension=2|3 got " << intrinsicDimension;

    if (!nef.is_simple())
      throw InvalidArgumentException(HERE) <<  "Nef polyhedron is not simple";

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
          // typename Polyhedron::Vertex_const_handle vch = vi;
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
          for (UnsignedInteger j = 0; j < dimension + 1; ++ j)
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
#if 0
            // filter out small simplex
            const Mesh simplexMesh(verticesI, IndicesCollection(Collection<Indices>(1, simplex)));
            if (!(simplexMesh.getVolume() > SpecFunc::Precision))
              continue;
#endif
            simplexColl.add(simplex);
          }
        }
        result.add(Mesh(verticesI, IndicesCollection(simplexColl)));
      }
    } // for nef.volumes
  } // 3d
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
}

/* Method load() reloads the object from the StorageManager */
void ConvexDecompositionMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
}


} /* namespace OTMESHING */
