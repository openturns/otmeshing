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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation.h>
#include <CGAL/Epick_d.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Triangulation2;
typedef CGAL::Delaunay_triangulation_3<K> Triangulation3;
typedef CGAL::Triangulation<CGAL::Epick_d< CGAL::Dynamic_dimension_tag > > TriangulationD;

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(CloudMesher);

static Factory<CloudMesher> Factory_CloudMesher;


/* Default constructor */
CloudMesher::CloudMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor method */
CloudMesher * CloudMesher::clone() const
{
  return new CloudMesher(*this);
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
    vertices = points;
    vertices.sortAccordingToAComponentInPlace(0);
    simplices = IndicesCollection(size - 1, dimension + 1);
    for (UnsignedInteger i = 0; i < size - 1; ++ i)
    {
      simplices(i, 0) = i;
      simplices(i, 1) = i + 1;
    }
  }
  else if (dimension == 2)
  {
    Triangulation2 triangulation;
    for (UnsignedInteger i = 0; i < size; ++ i)
      triangulation.insert({points(i, 0), points(i, 1)});

    UnsignedInteger iV = 0;
    std::map<Triangulation2::Point, UnsignedInteger> vertexToIndexMap; // cgal does not expose vertex indices
    for (Triangulation2::Vertex_handle vh : triangulation.finite_vertex_handles())
    {
      vertices.add(Point(vh->point().cartesian_begin(), vh->point().cartesian_end()));
      vertexToIndexMap[vh->point()] = iV;
      ++ iV;
    }

    simplices = IndicesCollection(triangulation.number_of_faces(), dimension + 1);
    UnsignedInteger iS = 0;
    for (Triangulation2::Face_handle f : triangulation.finite_face_handles())
    {
      auto triangle(triangulation.triangle(f));
      simplices(iS, 0) = vertexToIndexMap[triangle[0]];
      simplices(iS, 1) = vertexToIndexMap[triangle[1]];
      simplices(iS, 2) = vertexToIndexMap[triangle[2]];
      ++ iS;
    }
  }
  else if (dimension == 3)
  {
    Triangulation3 triangulation;
    for (UnsignedInteger i = 0; i < size; ++ i)
      triangulation.insert({points(i, 0), points(i, 1), points(i, 2)});

    UnsignedInteger iV = 0;
    std::map<Triangulation3::Point, UnsignedInteger> vertexToIndexMap; // cgal does not expose vertex indices
    for (Triangulation3::Vertex_handle vh : triangulation.finite_vertex_handles())
    {
      vertices.add(Point(vh->point().cartesian_begin(), vh->point().cartesian_end()));
      vertexToIndexMap[vh->point()] = iV;
      ++ iV;
    }

    simplices = IndicesCollection(triangulation.number_of_finite_cells(), dimension + 1);
    UnsignedInteger iS = 0;
    for (Triangulation3::Cell_handle f : triangulation.finite_cell_handles())
    {
      auto tetrahedron(triangulation.tetrahedron(f));
      simplices(iS, 0) = vertexToIndexMap[tetrahedron[0]];
      simplices(iS, 1) = vertexToIndexMap[tetrahedron[1]];
      simplices(iS, 2) = vertexToIndexMap[tetrahedron[2]];
      simplices(iS, 3) = vertexToIndexMap[tetrahedron[3]];
      ++ iS;
    }
  }
  else
  {
    TriangulationD triangulation(dimension);
    for (UnsignedInteger i = 0; i < size; ++ i)
    {
      const Point point(points[i]);
      triangulation.insert(TriangulationD::Point{dimension, point.begin(), point.end()});
    }

    UnsignedInteger iV = 0;
    std::map<TriangulationD::Point, UnsignedInteger> vertexToIndexMap; // cgal does not expose vertex indices
    for (TriangulationD::Vertex_iterator vi = triangulation.vertices_begin(); vi != triangulation.vertices_end(); ++ vi)
    {
      vertices.add(Point(vi->point().cartesian_begin(), vi->point().cartesian_end()));
      vertexToIndexMap[vi->point()] = iV;
      ++ iV;
    }

    simplices = IndicesCollection(triangulation.number_of_finite_full_cells(), dimension + 1);
    UnsignedInteger iS = 0;
    for (TriangulationD::Finite_full_cell_const_iterator cit = triangulation.finite_full_cells_begin(); cit != triangulation.finite_full_cells_end(); ++ cit)
    {
      for (UnsignedInteger j = 0; j < dimension + 1; ++ j)
      {
        TriangulationD::Vertex_handle vh = cit->vertex(j);
        simplices(iS, j) = vertexToIndexMap[vh->point()];
      }
      ++ iS;
    }
  }
  return Mesh(vertices, simplices);
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
  PersistentObject::save( adv );
}

/* Method load() reloads the object from the StorageManager */
void CloudMesher::load(Advocate & adv)
{
  PersistentObject::load( adv );
}


} /* namespace OTMESHING */
