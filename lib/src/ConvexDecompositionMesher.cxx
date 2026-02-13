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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

using namespace OT;

using KernelExact = CGAL::Exact_predicates_exact_constructions_kernel;
using Polyhedron = CGAL::Polyhedron_3<KernelExact>;
using HDS = Polyhedron::HalfedgeDS;
using Surface_mesh = CGAL::Surface_mesh<KernelExact::Point_3>;
using Nef_polyhedron = CGAL::Nef_polyhedron_3<KernelExact>;
using Point_3 = KernelExact::Point_3;

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
  if (dimension != 3)
    throw InvalidArgumentException(HERE) << "ConvexDecompositionMesher expected dimension=3 got " << dimension;

  const Sample vertices(mesh.getVertices());
  const IndicesCollection simplices(mesh.getSimplices());

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

    if (!nef.is_simple())
      throw InvalidArgumentException(HERE) <<  "Nef polyhedron is not simple";
  }
  else if (intrinsicDimension == 3)
  {
    const Point simplicesVolume(mesh.computeSimplicesVolume());

    // build from the volumetric mesh
    for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    {
      // LevelSetMesher can yield almost empty cells
      // possible workaround with key LevelSetMesher-SolveEquation=False
      if (!(simplicesVolume[i] > SpecFunc::Precision))
        continue;

      const UnsignedInteger i1 = simplices(i, 0);
      const UnsignedInteger i2 = simplices(i, 1);
      const UnsignedInteger i3 = simplices(i, 2);
      const UnsignedInteger i4 = simplices(i, 3);

      const Point_3 v1{vertices(i1, 0), vertices(i1, 1), vertices(i1, 2)};
      const Point_3 v2{vertices(i2, 0), vertices(i2, 1), vertices(i2, 2)};
      const Point_3 v3{vertices(i3, 0), vertices(i3, 1), vertices(i3, 2)};
      const Point_3 v4{vertices(i4, 0), vertices(i4, 1), vertices(i4, 2)};

      Polyhedron poly;
      poly.make_tetrahedron(v1, v2, v3, v4);

      const Nef_polyhedron tetra(poly);
      nef += tetra;
    }
  }
  else
    throw InvalidArgumentException(HERE) << "ConvexDecompositionMesher expected intrinsic dimension=2|3 got " << intrinsicDimension;

  CGAL::convex_decomposition_3(nef);

  Collection<Mesh> result;

  // Extract convex components

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
  }
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
