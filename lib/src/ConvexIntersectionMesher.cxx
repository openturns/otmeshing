//                                               -*- C++ -*-
/**
 *  @brief Intersection meshing algorithm for convex meshes
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
#include "openturns/SpecFunc.hxx"
#include "otmeshing/ConvexIntersectionMesher.hxx"
#include "otmeshing/CloudMesher.hxx"

#ifdef OPENTURNS_HAVE_CDDLIB
#include <setoper.h>
#include <cdd.h>
#endif

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(ConvexIntersectionMesher)
static const Factory<ConvexIntersectionMesher> Factory_ConvexIntersectionMesher;


/* Default constructor */
ConvexIntersectionMesher::ConvexIntersectionMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor */
ConvexIntersectionMesher * ConvexIntersectionMesher::clone() const
{
  return new ConvexIntersectionMesher(*this);
}

/* String converter */
String ConvexIntersectionMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << ConvexIntersectionMesher::GetClassName();
  return oss;
}

/* String converter */
String ConvexIntersectionMesher::__str__(const String & ) const
{
  return __repr__();
}

/* Here is the interface that all derived class must implement */

Mesh ConvexIntersectionMesher::build(const OT::Mesh & mesh1, const OT::Mesh & mesh2) const
{
  const UnsignedInteger dimension = mesh1.getDimension();
  if (mesh2.getDimension() != dimension)
    throw InvalidArgumentException(HERE) << "ConvexIntersectionMesher expected meshes of same dimension";

#ifdef OPENTURNS_HAVE_CDDLIB
  const Sample vertices1(mesh1.getVertices());
  const UnsignedInteger nv1 = vertices1.getSize();
  const Sample vertices2(mesh2.getVertices());
  const UnsignedInteger nv2 = vertices2.getSize();

  Sample vertices(0, dimension);
  CloudMesher cloudMesher;
  Collection<Indices> simplexColl;

  // initialize cddlib
  dd_ErrorType err = dd_NoError;
  dd_set_global_constants();

  // allocate V-representation
  dd_MatrixPtr m1 = dd_CreateMatrix(nv1, dimension + 1);
  dd_SetMatrixRepresentationType(m1, dd_Generator);
  dd_MatrixPtr m2 = dd_CreateMatrix(nv2, dimension + 1);
  dd_SetMatrixRepresentationType(m2, dd_Generator);
  for (UnsignedInteger j = 0; j < nv1; ++ j)
    // homogeneous coordinate
    dd_set_d(m1->matrix[j][0], 1.0);  
  for (UnsignedInteger j = 0; j < nv2; ++ j)
    // homogeneous coordinate
    dd_set_d(m2->matrix[j][0], 1.0);

  for (UnsignedInteger i1 = 0; i1 < nv1; ++ i1)
  {
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      dd_set_d(m1->matrix[i1][k + 1], vertices1(i1, k));
  }
  dd_PolyhedraPtr p1 = dd_DDMatrix2Poly(m1, &err);
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for mesh 1";

  // Convert V-representation to H-representation (inequalities)
  dd_MatrixPtr h1 = dd_CopyInequalities(p1);

  for (UnsignedInteger i2 = 0; i2 < nv2; ++ i2)
  {
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      dd_set_d(m2->matrix[i2][k + 1], vertices2(i2, k));
  } // Faces of the simplex
  dd_PolyhedraPtr p2 = dd_DDMatrix2Poly(m2, &err);
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for mesh 2";

  // Convert V-representation to H-representation (inequalities)
  dd_MatrixPtr h2 = dd_CopyInequalities(p2);
  dd_FreePolyhedra(p2);

  // Combine inequalities to compute intersection
  dd_MatrixAppendTo(&h2, h1);

  // Convert intersection back to V-representation
  dd_PolyhedraPtr intersectionV = dd_DDMatrix2Poly(h2, &err);
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for intersection";
  dd_FreeMatrix(h2);

  // retrieve vertices
  dd_MatrixPtr gen = dd_CopyGenerators(intersectionV);
  dd_FreePolyhedra(intersectionV);
  const UnsignedInteger verticesSize = vertices.getSize();
  const UnsignedInteger intersectionVerticesNumber = gen->rowsize; // empty intersection if zero
  if (intersectionVerticesNumber >= (dimension + 1))
  {
    // retrieve vertices
    Sample intersectionVertices(0, dimension);
    for (UnsignedInteger i = 0; i < intersectionVerticesNumber; ++i)
      {
        // First entry = 1 -> point, 0 -> ray
        if (dd_get_d(gen->matrix[i][0]) != 1.0)
          throw InternalException(HERE) << "assumed only points, no rays";

        Point vertex(dimension);
        for (UnsignedInteger j = 0; j < dimension; ++ j)
          vertex[j] = dd_get_d(gen->matrix[i][j + 1]);
        intersectionVertices.add(vertex);
      }

    // build simplices
    if (intersectionVerticesNumber == (dimension + 1))
    {
      // only one simplex
      Indices simplex(dimension + 1);
      simplex.fill(verticesSize);
      simplexColl.add(simplex);
      vertices.add(intersectionVertices);
    }
    else
    {
      // V>d+1, decompose into several simplices
      const Mesh intersectionMesh(cloudMesher.build(intersectionVertices));
      const IndicesCollection intersectionSimplices(intersectionMesh.getSimplices());
      const UnsignedInteger intersectionSimplicesNumber = intersectionMesh.getSimplicesNumber();
      for (UnsignedInteger i = 0; i < intersectionSimplicesNumber; ++ i)
      {
        Indices simplex(dimension + 1);
        for (UnsignedInteger j = 0; j <= dimension; ++ j)
          simplex[j] = verticesSize + intersectionSimplices(i, j);
        simplexColl.add(simplex);
      }
      // triangulation vertices are different than the ones in intersectionVertices
      vertices.add(intersectionMesh.getVertices());
    } // else V>d+1, decompose into several simplices
  } // if (intersectionVerticesNumber >= (dimension + 1))
  dd_FreeMatrix(h1);
  dd_FreePolyhedra(p1);

  // free cddlib objects
  dd_FreeMatrix(m1);
  dd_FreeMatrix(m2);
  dd_free_global_constants();
  
  // copy simplices
  IndicesCollection simplices(simplexColl.getSize(), dimension + 1);
  for (UnsignedInteger i = 0; i < simplexColl.getSize(); ++ i)
    for (UnsignedInteger j = 0; j <= dimension; ++ j)
      simplices(i, j) = simplexColl[i][j];

  return Mesh(vertices, simplices);
#else
  throw NotYetImplementedException(HERE) << "No cddlib support";
#endif
}

}
