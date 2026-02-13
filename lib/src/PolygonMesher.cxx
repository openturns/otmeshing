//                                               -*- C++ -*-
/**
 *  @brief Polygon meshing algorithm
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
#include "openturns/PersistentObjectFactory.hxx"
#include "otmeshing/PolygonMesher.hxx"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>


using KernelInexact = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = CGAL::Point_2<KernelInexact>;
using Polygon_2 = CGAL::Polygon_2<KernelInexact>;

struct Point2Hash
{
    std::size_t operator()(const Point_2& p) const
    {
        std::hash<double> h;
        return h(p.x()) ^ (h(p.y()) << 1);
    }
};

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(PolygonMesher)
static const Factory<PolygonMesher> Factory_PolygonMesher;


/* Default constructor */
PolygonMesher::PolygonMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor */
PolygonMesher * PolygonMesher::clone() const
{
  return new PolygonMesher(*this);
}

/* String converter */
String PolygonMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << PolygonMesher::GetClassName();
  return oss;
}

/* String converter */
String PolygonMesher::__str__(const String & ) const
{
  return __repr__();
}

/* Here is the interface that all derived class must implement */

Mesh PolygonMesher::build(const Sample & points) const
{
  const UnsignedInteger dimension = points.getDimension();
  const UnsignedInteger size = points.getSize();

  if (size < 3)
    throw InvalidArgumentException(HERE) << "PolygonMesher expected points of size >=3, got " << size;

  const Point stddev(points.computeStandardDeviation());
  Indices intrinsic;
  for (UnsignedInteger j = 0; j < dimension; ++ j)
    if (stddev[j] > 0.0)
      intrinsic.add(j);
  if (intrinsic.getSize() != 2)
    throw InvalidArgumentException(HERE) << "PolygonMesher expected an intrinsic dimension of 2, got " << size;

  // the indices are not stored so we need to query a map
  std::vector<Point_2> points2;
  std::unordered_map<Point_2, UnsignedInteger> vertexToIndexMap;
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    const Point_2 p2(points(i, intrinsic[0]), points(i, intrinsic[1]));
    points2.push_back(p2);
    // A polygon is called simple if there is no pair of nonconsecutive edges sharing a point
    if (vertexToIndexMap.count(p2) > 0)
      throw InvalidArgumentException(HERE) << "PolygonMesher expected a simple polygon (no redundant vertex)";
    vertexToIndexMap[p2] = i;
  }

  const Polygon_2 polygon(points2.begin(), points2.end());

  // perform triangulation
  std::vector<Polygon_2> triangles;
  const CGAL::Polygon_triangulation_decomposition_2<KernelInexact> decomposition2;
  decomposition2(polygon, std::back_inserter(triangles));

  // retrieve triangles
  IndicesCollection simplices(triangles.size(), dimension + 1);
  for (UnsignedInteger i = 0; i < triangles.size(); ++ i)
  {
    for (UnsignedInteger j = 0; j < 3; ++ j)
    {
      const Point_2 & vh = triangles[i].vertex(j);
      simplices(i, j) = vertexToIndexMap[vh];
    }
    // repeat last index to mark intrinsic dimension
    for (UnsignedInteger j = 3; j <= dimension; ++ j)
      simplices(i, j) = simplices(i, 2);
  }
  return Mesh(points, simplices);
}

/* Method save() stores the object through the StorageManager */
void PolygonMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
}

/* Method load() reloads the object from the StorageManager */
void PolygonMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
}

}
