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
#include "otmeshing/CloudMesher.hxx"
#include <openturns/PersistentObjectFactory.hxx>

#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/Delaunay_triangulation.h>

#if 0
#include "libqhull_r/qhull_ra.h"
#endif

using DefaultTriangulation = CGAL::Triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag> > ;
using DelaunayTriangulation = CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag> >;

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(CloudMesher)

static Factory<CloudMesher> Factory_CloudMesher;


/* Default constructor */
CloudMesher::CloudMesher(const TriangulationMethod method)
  : PersistentObject()
  , triangulationMethod_(method)
{
  // Nothing to do
}

/* Virtual constructor method */
CloudMesher * CloudMesher::clone() const
{
  return new CloudMesher(*this);
}


template <class TriangulationType>
Mesh buildTriangulation(const Sample & points)
{
  const UnsignedInteger dimension = points.getDimension();
  const UnsignedInteger size = points.getSize();
  TriangulationType triangulation(dimension);
  std::vector<typename TriangulationType::Point> pts(size);
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    const Point point(points[i]);
    pts[i] = typename TriangulationType::Point{dimension, point.begin(), point.end()};
  }
  // it is much faster to insert vertices by batch
  triangulation.insert(pts.begin(), pts.end());

  // the vertices are reordered by the triangulation
  std::unordered_map<typename TriangulationType::Vertex_iterator, UnsignedInteger> vertexToIndexMap;
  UnsignedInteger vertexIndex = 0;
  Sample vertices(0, dimension);

  // the infinite first vertex can be skipped
  for (typename TriangulationType::Vertex_iterator vi = ++ triangulation.vertices_begin(); vi != triangulation.vertices_end(); ++ vi)
  {
    vertices.add(Point(vi->point().cartesian_begin(), vi->point().cartesian_end()));
    vertexToIndexMap[vi] = vertexIndex;
    ++ vertexIndex;
  }

  IndicesCollection simplices(triangulation.number_of_finite_full_cells(), dimension + 1);
  UnsignedInteger simplexIndex = 0;
  for (typename TriangulationType::Finite_full_cell_const_iterator cit = triangulation.finite_full_cells_begin(); cit != triangulation.finite_full_cells_end(); ++ cit)
  {
    for (UnsignedInteger j = 0; j < dimension + 1; ++ j)
    {
      const typename TriangulationType::Vertex_handle vh = cit->vertex(j);
      simplices(simplexIndex, j) = vertexToIndexMap[vh];
    }
    ++ simplexIndex;
  }
  return Mesh(vertices, simplices);
}


Mesh CloudMesher::build(const Sample & points) const
{
  const UnsignedInteger dimension = points.getDimension();
  const UnsignedInteger size = points.getSize();
  if (!dimension)
    throw InvalidArgumentException(HERE) << "CloudMesher expected a non-null dimension";
  if (size < dimension + 1)
    throw InvalidArgumentException(HERE) << "CloudMesher expected a size of at least " << dimension + 1 << " got " << size;
  Sample vertices(0, points.getDimension());
  if (dimension == 1)
  {
    // special case for dim=1 to avoid special handling in the generic part
    vertices.add(points.getMin());
    vertices.add(points.getMax());
    IndicesCollection simplices(1, dimension + 1);
    simplices(0, 1) = 1;
    return Mesh(vertices, simplices);
  }

#if 0
  QHULL_LIB_CHECK

  // Create qhull context
  qhT qh_qh;
  qhT *qh = &qh_qh;
  qh_zero(qh, stderr);

  // Run Qhull
  const String qhull_cmd("qhull d Qt Qx Qz"); // options: delaunay + triangulated output + deterministic output + infinity point
  int rc = qh_new_qhull(qh, dimension, size,
                        const_cast<Scalar*>(points.getImplementation()->data()),
                        False, /* ismalloc */
                        const_cast<char*>(qhull_cmd.c_str()),
                        NULL, NULL);

  if (rc != 0)
  {
    qh_freeqhull(qh, !qh_ALL);
    throw InternalException(HERE) << "qh_new_qhull exit code: " << rc;
  }

  // build the vertices
  Indices inputIndexToHullIndex(size, size);
  vertexT *vertex = NULL, **vertexp = NULL;
  UnsignedInteger i = 0;
  FORALLvertices
  {
    if (!vertex->deleted)
    {
      // qh_pointid gives indices wrt the original input sample
      const SignedInteger inputIdx = qh_pointid(qh, vertex->point);

      // infinite vertex (Qz option)
      if (static_cast<UnsignedInteger>(inputIdx) == size)
        continue;

      Point point(dimension);
      // assume vertex->point is an array of double
      std::copy(vertex->point, vertex->point + dimension, point.begin());
      vertices.add(point);

      inputIndexToHullIndex[inputIdx] = i;
      ++ i;
    }
  }

  // build the faces
  Collection<Indices> simplexColl;
  Indices used(vertices.getSize());
  facetT *facet = NULL;
  FORALLfacets
  {
    if (!facet->upperdelaunay) /* skip upper facets in 2D */
    {
      Indices simplex(dimension + 1);
      UnsignedInteger j = 0;
      FOREACHvertex_(facet->vertices)
      {
        const SignedInteger hullIdx = qh_pointid(qh, vertex->point);
        simplex[j] = inputIndexToHullIndex[hullIdx];
        used[simplex[j]] = 1;
        ++ j;
      }

      simplexColl.add(simplex);
    }
  }

  // cleanup
  qh_freeqhull(qh, !qh_ALL);
  int curlong = 0, totlong = 0;
  qh_memfreeshort(qh, &curlong, &totlong);
  if (curlong || totlong)
    throw InternalException(HERE) << "qh_memfreeshort: did not free " << totlong <<" bytes (" << curlong << " blocks)";

  return Mesh(vertices, IndicesCollection(simplexColl));
#else
  switch (triangulationMethod_)
  {
    case BASIC:
      return buildTriangulation<DefaultTriangulation>(points);
    case DELAUNAY:
      return buildTriangulation<DelaunayTriangulation>(points);
    default:
      throw InvalidArgumentException(HERE) << "Unknown triangulation method: " << triangulationMethod_;
  }
#endif
}

/* String converter */
String CloudMesher::__repr__() const
{
  OSS oss;
  oss << "class=" << CloudMesher::GetClassName();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void CloudMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("triangulationMethod_", triangulationMethod_);
}

/* Method load() reloads the object from the StorageManager */
void CloudMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("triangulationMethod_", triangulationMethod_);
}


} /* namespace OTMESHING */
