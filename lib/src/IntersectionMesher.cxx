//                                               -*- C++ -*-
/**
 *  @brief Intersection meshing algorithm
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
#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/SpecFunc.hxx>

#include "otmeshing/IntersectionMesher.hxx"
#include "otmeshing/CloudMesher.hxx"
#include "otmeshing/ConvexDecompositionMesher.hxx"
#include "otmeshing/UnionMesher.hxx"

#ifdef OPENTURNS_HAVE_CDDLIB
#include <setoper.h>
#include <cdd.h>
#endif

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(IntersectionMesher)
static const Factory<IntersectionMesher> Factory_IntersectionMesher;


/* Default constructor */
IntersectionMesher::IntersectionMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor */
IntersectionMesher * IntersectionMesher::clone() const
{
  return new IntersectionMesher(*this);
}

/* String converter */
String IntersectionMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << IntersectionMesher::GetClassName();
  return oss;
}

/* String converter */
String IntersectionMesher::__str__(const String & ) const
{
  return __repr__();
}

Mesh IntersectionMesher::build(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

  const UnsignedInteger dimension = coll[0].getDimension();
  if ((coll.getSize() == 2) && (dimension == 3000000)) // TODO: enable this
  {
    ConvexDecompositionMesher convexDecompositionMesher;
    const Collection<Mesh> decomposition1(convexDecompositionMesher.build(coll[0]));
    const Collection<Mesh> decomposition2(convexDecompositionMesher.build(coll[1]));
    Sample vertices(0, dimension);
    Collection<Indices> simplexColl;
    for (UnsignedInteger i1 = 0; i1 < decomposition1.getSize(); ++ i1)
    {
      for (UnsignedInteger i2 = 0; i2 < decomposition2.getSize(); ++ i2)
      {
        // TODO: parallelize ?
        const Mesh intersection12(buildConvex(Collection<Mesh>({decomposition1[i1], decomposition2[i2]})));
        if (intersection12.getSimplicesNumber() == 0)
          continue;

        vertices.add(intersection12.getVertices());
        const IndicesCollection simplicesI(intersection12.getSimplices());
        for (UnsignedInteger k = 0; k < simplicesI.getSize(); ++ k)
        {
          simplexColl.add(Indices(simplicesI.getImplementation()->cbegin_at(k),
                                  simplicesI.getImplementation()->cbegin_at(k) + dimension + 1));
        }
      }
    }
    Mesh result(vertices, IndicesCollection(simplexColl));
    if (recompress_)
      result = UnionMesher::CompressMesh(result);
    return result;
  } // dim=3

  Collection<Mesh> todo(coll);
  while (todo.getSize() > 1)
  {
    Collection<Mesh> done(todo.getSize() / 2);
    // TODO: parallelize ?
    for (UnsignedInteger i = 0; i < todo.getSize() / 2; ++ i)
      done[i] = build2(todo[2 * i], todo[2 * i + 1]);

    // report odd element
    if (todo.getSize() % 2)
      done.add(todo[todo.getSize() - 1]);

    todo = done;
  }
  return todo[0];
}

#ifdef OPENTURNS_HAVE_CDDLIB
String cdd_error_to_string(const dd_ErrorType err)
{
  switch (err)
  {
    case dd_DimensionTooLarge:
      return "Dimension too large";
    case dd_ImproperInputFormat:
      return "Improper input format";
    case dd_NegativeMatrixSize:
      return "Negative matrix size";
    case dd_EmptyVrepresentation:
      return "Empty V-representation";
    case dd_EmptyHrepresentation:
      return "Empty H-representation";
    case dd_EmptyRepresentation:
      return "Empty representation";
    case dd_IFileNotFound:
      return "Input file not found";
    case dd_OFileNotOpen:
      return "Output file not open";
    case dd_NoLPObjective:
      return "No LP objective specified";
    case dd_NoRealNumberSupport:
      return "No real number support (library built without GMP?)";
    case dd_NotAvailForH:
      return "Operation not available for H-representation";
    case dd_NotAvailForV:
      return "Operation not available for V-representation";
    case dd_CannotHandleLinearity:
      return "Cannot handle linearity in this context";
    case dd_RowIndexOutOfRange:
      return "Row index out of range";
    case dd_ColIndexOutOfRange:
      return "Column index out of range";
    case dd_LPCycling:
      return "LP cycling detected";
    case dd_NumericallyInconsistent:
      return "Numerical inconsistency detected";
    case dd_NoError:
      return "No error";
    default:
        return "Unknown cddlib error";
  }
}
#endif

Mesh IntersectionMesher::build2(const Mesh & mesh1, const Mesh & mesh2) const
{
  const UnsignedInteger dimension = mesh1.getDimension();
  if (mesh2.getDimension() != dimension)
    throw InvalidArgumentException(HERE) << "IntersectionMesher expected meshes of same dimension";

#ifdef OPENTURNS_HAVE_CDDLIB
  const IndicesCollection simplices1(mesh1.getSimplices());
  const IndicesCollection simplices2(mesh2.getSimplices());
  const Sample vertices1(mesh1.getVertices());
  const Sample vertices2(mesh2.getVertices());
  const UnsignedInteger ns1 = mesh1.getSimplicesNumber();
  const UnsignedInteger ns2 = mesh2.getSimplicesNumber();

  CloudMesher cloudMesher;
  Collection<Mesh> intersectionColl;

  // initialize cddlib
  dd_ErrorType err = dd_NoError;
  dd_set_global_constants();

  // allocate V-representation
  dd_MatrixPtr m1 = dd_CreateMatrix(dimension + 1, dimension + 1);
  dd_SetMatrixRepresentationType(m1, dd_Generator);
  dd_MatrixPtr m2 = dd_CreateMatrix(dimension + 1, dimension + 1);
  dd_SetMatrixRepresentationType(m2, dd_Generator);
  for (UnsignedInteger j = 0; j <= dimension; ++ j)
  {
    // homogeneous coordinate
    dd_set_d(m1->matrix[j][0], 1.0);  
    dd_set_d(m2->matrix[j][0], 1.0);
  }

  // mesh1 simplices loop
  for (UnsignedInteger i1 = 0; i1 < ns1; ++ i1)
  {
    Point lower1(dimension, SpecFunc::Infinity);
    Point upper1(dimension, -SpecFunc::Infinity);
    // build V-representation of simplex
    for (UnsignedInteger j = 0; j <= dimension; ++ j)
    {
      const UnsignedInteger vi1j = simplices1(i1, j);
      for (UnsignedInteger k = 0; k < dimension; ++ k)
      {
        dd_set_d(m1->matrix[j][k + 1], vertices1(vi1j, k));
        lower1[k] = std::min(lower1[k], vertices1(vi1j, k));
        upper1[k] = std::max(upper1[k], vertices1(vi1j, k));
      }
    }
    dd_PolyhedraPtr p1 = dd_DDMatrix2Poly(m1, &err);
    if (err != dd_NoError)
      throw InternalException(HERE) << "dd_DDMatrix2Poly failed i1=" << i1 << ": " << cdd_error_to_string(err);

    // Convert V-representation to H-representation (inequalities)
    dd_MatrixPtr h1 = dd_CopyInequalities(p1);

    // mesh2 simplices loop
    for (UnsignedInteger i2 = 0; i2 < ns2; ++ i2)
    {
      Point lower2(dimension, SpecFunc::Infinity);
      Point upper2(dimension, -SpecFunc::Infinity);
      // build V-representation of simplex
      for (UnsignedInteger j = 0; j <= dimension; ++ j)
      {
        const UnsignedInteger vi2j = simplices2(i2, j);
        for (UnsignedInteger k = 0; k < dimension; ++ k)
        {
          dd_set_d(m2->matrix[j][k + 1], vertices2(vi2j, k));
          lower2[k] = std::min(lower2[k], vertices2(vi2j, k));
          upper2[k] = std::max(upper2[k], vertices2(vi2j, k));
        }
      }
      Bool toSkip = false;
      for (UnsignedInteger k = 0; k < dimension; ++ k)
      {
        toSkip = std::max(lower1[k], lower2[k]) >= std::min(upper1[k], upper2[k]);
        if (toSkip)
          break;
      }
      if (toSkip)
        continue;

      dd_PolyhedraPtr p2 = dd_DDMatrix2Poly(m2, &err);
      if (err != dd_NoError)
        throw InternalException(HERE) << "dd_DDMatrix2Poly failed i2=" << i2 << ": " << cdd_error_to_string(err);

      // Convert V-representation to H-representation (inequalities)
      dd_MatrixPtr h2 = dd_CopyInequalities(p2);
      dd_FreePolyhedra(p2);

      // Combine inequalities to compute intersection
      dd_MatrixAppendTo(&h2, h1);
      dd_SetMatrixRepresentationType(h2, dd_Inequality);

      // Convert intersection back to V-representation
      dd_PolyhedraPtr intersectionV = dd_DDMatrix2Poly(h2, &err);
      if (err != dd_NoError)
        throw InternalException(HERE) << "dd_DDMatrix2Poly failed for intersection: "  << cdd_error_to_string(err);
      dd_FreeMatrix(h2);

      // retrieve vertices
      dd_MatrixPtr gen = dd_CopyGenerators(intersectionV);
      dd_FreePolyhedra(intersectionV);
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
          simplex.fill(); // orientation may be incorrect
          const IndicesCollection intersectionSimplices(Collection<Indices>(1, simplex));
          intersectionColl.add(Mesh(intersectionVertices, intersectionSimplices));
        }
        else
        {
          // V>d+1, decompose into several simplices
          const Mesh intersectionMesh(cloudMesher.build(intersectionVertices));
          intersectionColl.add(intersectionMesh);
        }
      } // if (intersectionVerticesNumber >= (dimension + 1))
    } // mesh2 simplices loop
    dd_FreeMatrix(h1);
    dd_FreePolyhedra(p1);

  } // mesh1 simplices loop

  // free cddlib objects
  dd_FreeMatrix(m1);
  dd_FreeMatrix(m2);
  dd_free_global_constants();

  Mesh result(UnionMesher().build(intersectionColl));
  if (recompress_)
    result = UnionMesher::CompressMesh(result);
  return result;
#else
  throw NotYetImplementedException(HERE) << "No cddlib support";
#endif
}


Mesh IntersectionMesher::buildConvex(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

  const UnsignedInteger dimension = coll[0].getDimension();
  for (UnsignedInteger i = 1; i < size; ++ i)
    if (coll[i].getDimension() != dimension)
      throw InvalidArgumentException(HERE) << "IntersectionMesher expected meshes of same dimension";

#ifdef OPENTURNS_HAVE_CDDLIB
  CloudMesher cloudMesher;
  Collection<Mesh> intersectionColl;
  Point lower1(dimension, -SpecFunc::Infinity);
  Point upper1(dimension, SpecFunc::Infinity);
  UnsignedInteger prunedNumber = 0;

  // initialize cddlib
  dd_ErrorType err = dd_NoError;
  dd_set_global_constants();

  // allocate H-representation of intersection
  dd_MatrixPtr intersectionH = dd_CreateMatrix(0, dimension + 1);
  dd_SetMatrixRepresentationType(intersectionH, dd_Inequality);

  // for each convex
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    const Sample vertices1(coll[i].getVertices());
    const UnsignedInteger nv1 = vertices1.getSize();

    // bbox pruning
    const Point min1(vertices1.getMin());
    const Point max1(vertices1.getMax());
    Bool toSkip = false;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      toSkip = std::max(lower1[k], min1[k]) >= std::min(upper1[k], max1[k]);
      if (toSkip)
        break;
    }
    if (toSkip)
    {
      ++ prunedNumber;
      continue;
    }
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      lower1[k] = std::max(lower1[k], min1[k]);
      upper1[k] = std::min(upper1[k], max1[k]);
    }

    // allocate V-representation
    dd_MatrixPtr m1 = dd_CreateMatrix(nv1, dimension + 1);
    dd_SetMatrixRepresentationType(m1, dd_Generator);
    for (UnsignedInteger i1 = 0; i1 < nv1; ++ i1)
    {
      // homogeneous coordinate
      dd_set_d(m1->matrix[i1][0], 1.0);
      for (UnsignedInteger k = 0; k < dimension; ++ k)
      {
        dd_set_d(m1->matrix[i1][k + 1], vertices1(i1, k));
      }
    }

    dd_PolyhedraPtr p1 = dd_DDMatrix2Poly(m1, &err);
    if (err != dd_NoError)
      throw InternalException(HERE) << "dd_DDMatrix2Poly failed for mesh 1: " << cdd_error_to_string(err);

    // Convert V-representation to H-representation (inequalities)
    dd_MatrixPtr h1 = dd_CopyInequalities(p1);

    // Combine inequalities
    dd_MatrixAppendTo(&intersectionH, h1);

    // free memory
    dd_FreeMatrix(m1);
    dd_FreePolyhedra(p1);
    dd_FreeMatrix(h1);

  } // i loop

  // empty intersection
  if (size - prunedNumber == 1)
    return Mesh(Sample(0, dimension));

  // Convert intersection back to V-representation
  dd_PolyhedraPtr intersectionV = dd_DDMatrix2Poly(intersectionH, &err);
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for intersection: " << cdd_error_to_string(err);
  dd_FreeMatrix(intersectionH);

  // retrieve vertices
  dd_MatrixPtr gen = dd_CopyGenerators(intersectionV);
  dd_FreePolyhedra(intersectionV);
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
      simplex.fill(); // orientation may be incorrect
      const IndicesCollection intersectionSimplices(Collection<Indices>(1, simplex));
      intersectionColl.add(Mesh(intersectionVertices, intersectionSimplices));
    }
    else
    {
      // V>d+1, decompose into several simplices
      const Mesh intersectionMesh(cloudMesher.build(intersectionVertices));
      intersectionColl.add(intersectionMesh);
    }
  } // if (intersectionVerticesNumber >= (dimension + 1))

  dd_FreeMatrix(gen);
  dd_free_global_constants();

  const Mesh result(UnionMesher().build(intersectionColl));
  return result;
#else
  throw NotYetImplementedException(HERE) << "No cddlib support";
#endif
}

Mesh IntersectionMesher::buildCylinder(const Collection<Cylinder> & coll) const
{
  CloudMesher mesher;
  const UnsignedInteger size = coll.getSize();
  Collection<Mesh> collMesh(size);
  for (UnsignedInteger i = 0; i < size; ++ i)
    collMesh[i] = mesher.build(coll[i].getVertices());
  return buildConvex(collMesh);
}

/* Recompression flag accessor */
void IntersectionMesher::setRecompress(const Bool recompress)
{
  recompress_ = recompress;
}

Bool IntersectionMesher::getRecompress() const
{
  return recompress_;
}

/* Method save() stores the object through the StorageManager */
void IntersectionMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("recompress_", recompress_);
}

/* Method load() reloads the object from the StorageManager */
void IntersectionMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("recompress_", recompress_);
}

}
