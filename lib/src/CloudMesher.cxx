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
  IndicesCollection simplices;
  if (dimension == 1)
  {
    // special case for dim=1 to avoid special handling in the generic part
    vertices.add(points.getMin());
    vertices.add(points.getMax());
    simplices = IndicesCollection(1, dimension + 1);
    simplices(0, 1) = 1;
    return Mesh(vertices, simplices);
  }
  else
  {
    switch (triangulationMethod_)
    {
      case BASIC:
        return buildTriangulation<DefaultTriangulation>(points);
      case DELAUNAY:
        return buildTriangulation<DelaunayTriangulation>(points);
      default:
        throw InvalidArgumentException(HERE) << "Unknown triangulation method: " << triangulationMethod_;
    }
  }
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
