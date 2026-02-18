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
#include "otmeshing/ConvexHullMesher.hxx"
#include <openturns/PersistentObjectFactory.hxx>

#include <CGAL/Epick_d.h>
#include <CGAL/Triangulation.h>

#ifdef OPENTURNS_HAVE_QHULL
#include "libqhull_r/qhull_ra.h"
#endif

using DefaultTriangulation = CGAL::Triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag> >;

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(ConvexHullMesher)

static Factory<ConvexHullMesher> Factory_ConvexHullMesher;


/* Default constructor */
ConvexHullMesher::ConvexHullMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor method */
ConvexHullMesher * ConvexHullMesher::clone() const
{
  return new ConvexHullMesher(*this);
}


template <class TriangulationType>
Mesh buildConvexHull(const Sample & points)
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

  // retrieve the number of full cells incident to the vertex at infinity (ie the number of infinite cells)
  UnsignedInteger facetNumber = 0;
  for (typename TriangulationType::Full_cell_iterator cit = triangulation.full_cells_begin(); cit != triangulation.full_cells_end(); ++ cit)
  {
    if (triangulation.is_infinite(cit))
      ++ facetNumber;
  }

  IndicesCollection simplices(facetNumber, dimension + 1);
  UnsignedInteger facetIndex = 0;
  Sample vertices(0, dimension);
  UnsignedInteger vertexIndex = 0;
  for (typename TriangulationType::Full_cell_iterator cit = triangulation.full_cells_begin(); cit != triangulation.full_cells_end(); ++ cit)
  {
    if (triangulation.is_infinite(cit))
    {
      // retrieve the facet incident to the infinite vertex
      const typename TriangulationType::Facet ft(cit, cit->index(triangulation.infinite_vertex()));

      // retrieve the cell vertex index outside the facet
      const UnsignedInteger covertexIndex = triangulation.index_of_covertex(ft);

      // this index skips covertexIndex and covers [0; dimension-1]
      UnsignedInteger j2 = 0;

      // for each vertex index in the cell
      for (UnsignedInteger j = 0; j < dimension + 1; ++ j)
      {
        if (j != covertexIndex)
        {
          const typename TriangulationType::Vertex_handle vh = cit->vertex(j);

          // check if this vertex was visited yet
          if (!vertexToIndexMap.count(vh))
          {
            vertices.add(Point(vh->point().cartesian_begin(), vh->point().cartesian_end()));
            vertexToIndexMap[vh] = vertexIndex;
            ++ vertexIndex;
          }

          simplices(facetIndex, j2) = vertexToIndexMap[vh];
          ++ j2;
        }
      }
      // repeat the last index to set the intrinsic dimension to be dimension - 1
      simplices(facetIndex, dimension) = simplices(facetIndex, dimension - 1);

      ++ facetIndex;
    }
  }

  return Mesh(vertices, simplices);
}


Mesh ConvexHullMesher::build(const Sample & points) const
{
  const UnsignedInteger dimension = points.getDimension();
  const UnsignedInteger size = points.getSize();
  if (!dimension)
    throw InvalidArgumentException(HERE) << "ConvexHullMesher expected a non-null dimension";
  if (size < dimension + 1)
    throw InvalidArgumentException(HERE) << "ConvexHullMesher expected a size of at least " << dimension + 1 << " got " << size;
  Sample vertices(0, dimension);
  if (dimension == 1)
  {
    // special case for dim=1 to avoid special handling in the generic part
    vertices.add(points.getMin());
    vertices.add(points.getMax());
    IndicesCollection simplices(1, dimension + 1);
    simplices(0, 1) = 1;
    return Mesh(vertices, simplices);
  }
#ifdef OPENTURNS_HAVE_QHULL
  QHULL_LIB_CHECK

  // Create qhull context
  qhT qh_qh;
  qhT *qh = &qh_qh;
  qh_zero(qh, stderr);

  // Run Qhull
  const String qhull_cmd("qhull Qt Qx"); // options: triangulated output + deterministic output
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
        ++ j;
      }

      // repeat the last value to mark the intrinsic dimension
      simplex[dimension] = simplex[dimension - 1];
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
  return buildConvexHull<DefaultTriangulation>(points);
#endif
}

/* String converter */
String ConvexHullMesher::__repr__() const
{
  OSS oss;
  oss << "class=" << ConvexHullMesher::GetClassName();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void ConvexHullMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
}

/* Method load() reloads the object from the StorageManager */
void ConvexHullMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
}


} /* namespace OTMESHING */
