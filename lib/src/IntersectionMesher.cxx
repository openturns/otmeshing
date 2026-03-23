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
#include <chrono>

#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/SpecFunc.hxx>
#include <openturns/TBBImplementation.hxx>

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

struct IntersectionMesherConvexSamplePolicy
{
  const IntersectionMesher & intersectionMesher_;
  const Collection<Sample> & input1_;
  const Collection<Sample> & input2_;
  UnsignedInteger done_;
  Collection<Sample> & output_;
  UnsignedInteger stride_;

  IntersectionMesherConvexSamplePolicy(const IntersectionMesher & intersectionMesher,
                                      const Collection<Sample> & input1,
                                      const Collection<Sample> & input2,
                                      const UnsignedInteger done,
                                      Collection<Sample> & output)
    : intersectionMesher_(intersectionMesher)
    , input1_(input1)
    , input2_(input2)
    , done_(done)
    , output_(output)
    , stride_(input1.getSize())
  {}

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger n = r.begin(); n != r.end(); ++n)
    {
      const UnsignedInteger i = (n + done_) % stride_;
      const UnsignedInteger j = (n + done_) / stride_;
      output_[n] = intersectionMesher_.buildConvexSample({input1_[i], input2_[j]});
    }
  }
};

Mesh IntersectionMesher::build(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

  // build decomposition of first mesh
  ConvexDecompositionMesher convexDecompositionMesher;
  convexDecompositionMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
#if 0
  LOGTRACE(OSS() << "Build decomposition of mesh 0");
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
#endif
  Collection<Sample> unionCurrent;
  const Collection<Mesh> baseDecomposition0(convexDecompositionMesher.build(coll[0]));
  const UnsignedInteger baseDecompositionSize0 = baseDecomposition0.getSize();
  for (UnsignedInteger j = 0; j < baseDecompositionSize0; ++ j)
    unionCurrent.add(baseDecomposition0[j].getVertices());
#if 0
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  const Scalar timeDuration = std::chrono::duration<Scalar>(t1 - t0).count();
  LOGTRACE(OSS() << "Got " << baseDecompositionSize0 << " parts in t=" << timeDuration << "s");
#endif

  // for each remaining mesh i
  for (UnsignedInteger i = 1; i < size; ++i)
  {
    // build decomposition
#if 0
    LOGTRACE(OSS() << "Build decomposition of mesh " << i);
    t0 = std::chrono::steady_clock::now();
#endif
    Collection<Sample> unionNext;
    const Collection<Mesh> baseDecomposition(convexDecompositionMesher.build(coll[i]));
    const UnsignedInteger baseDecompositionSize = baseDecomposition.getSize();
    //LOGTRACE(OSS() << "Got " << baseDecompositionSize << " parts");
    for (UnsignedInteger j = 0; j < baseDecompositionSize; ++ j)
      unionNext.add(baseDecomposition[j].getVertices());
#if 0
    t1 = std::chrono::steady_clock::now();
    timeDuration = std::chrono::duration<Scalar>(t1 - t0).count();
    LOGTRACE(OSS() << "Got " << baseDecompositionSize << " parts in t=" << timeDuration << "s");
    t0 = std::chrono::steady_clock::now();
#endif
    // loop over intersections
    const UnsignedInteger toDoSize = unionCurrent.getSize() * unionNext.getSize();
    const UnsignedInteger blockSize = ResourceMap::GetAsUnsignedInteger("IntersectionMesher-BlockSize");
    // LOGTRACE(OSS() << "Ready to compute " << toDoSize << " pairwise intersections by blocks of size " << blockSize);
    Collection<Sample> result(0);
    Collection<Sample> resultChunk(blockSize);
    for (UnsignedInteger done = 0; done < toDoSize; done += blockSize)
    {
      const UnsignedInteger actualBlockSize = std::min(blockSize, toDoSize - done);
      const IntersectionMesherConvexSamplePolicy policy(*this, unionCurrent, unionNext, done, resultChunk);
      TBBImplementation::ParallelFor(0, actualBlockSize, policy);

      // prune empty intersections
      for (UnsignedInteger i0 = 0; i0 < actualBlockSize; ++ i0)
        if (resultChunk[i0].getSize())
          result.add(resultChunk[i0]);
    }
#if 0
    t1 = std::chrono::steady_clock::now();
    timeDuration = std::chrono::duration<Scalar>(t1 - t0).count();
    LOGTRACE(OSS() << "Done, t=" << timeDuration << "s");
#endif
    // early exit if there are no non-empty intersections at this stage
    unionCurrent = result;
    if (!unionCurrent.getSize())
      return Mesh(Sample(0, coll[i].getDimension()));
  } // for mesh i

  // build mesh of union of remaining intersections
  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(unionCurrent.getSize());
  for (UnsignedInteger i = 0; i < unionCurrent.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(unionCurrent[i]);
  return UnionMesher().build(collMesh);
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


Mesh IntersectionMesher::buildConvex(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];
  Collection<Sample> collS(size);
  for (UnsignedInteger i = 0; i < size; ++ i)
    collS[i] = coll[i].getVertices();
  const Sample intersectionVertices(buildConvexSample(collS));
  const UnsignedInteger dimension = coll[0].getDimension();
  const UnsignedInteger intersectionVerticesNumber = intersectionVertices.getSize();

  // build mesh
  Mesh result;
  if (intersectionVerticesNumber == (dimension + 1))
  {
    // only one simplex
    Indices simplex(dimension + 1);
    simplex.fill(); // orientation may be incorrect
    const IndicesCollection intersectionSimplices(Collection<Indices>(1, simplex));
    result = Mesh(intersectionVertices, intersectionSimplices);
  }
  else if (intersectionVerticesNumber > (dimension + 1))
  {
    // V>d+1, decompose into several simplices
    result = CloudMesher().build(intersectionVertices);
  }
  return result;
}


Sample IntersectionMesher::buildConvexSample(const Collection<Sample> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Sample(0, 0);
  else if (size == 1)
    return coll[0];

  const UnsignedInteger dimension = coll[0].getDimension();
  for (UnsignedInteger i = 1; i < size; ++ i)
    if (coll[i].getDimension() != dimension)
      throw InvalidArgumentException(HERE) << "IntersectionMesher expected vertices of same dimension";

  Sample result(0, dimension);
  Point lower1(dimension, -SpecFunc::Infinity);
  Point upper1(dimension, SpecFunc::Infinity);

  // bbox pruning
  Indices remainingIndices;
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    const Sample vertices1(coll[i]);
    const Point min1(vertices1.getMin());
    const Point max1(vertices1.getMax());
    Bool pruned = false;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      pruned = std::max(lower1[k], min1[k]) >= std::min(upper1[k], max1[k]);
      if (pruned)
        break;
    }
    if (pruned)
      continue;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      lower1[k] = std::max(lower1[k], min1[k]);
      upper1[k] = std::min(upper1[k], max1[k]);
    }
    remainingIndices.add(i);
  }

  const UnsignedInteger remainingSize = remainingIndices.getSize();
  if (remainingSize == 1)
    return result;

#ifdef OPENTURNS_HAVE_CDDLIB

  // initialize cddlib
  dd_ErrorType err = dd_NoError;

  // allocate H-representation of intersection
  dd_MatrixPtr intersectionH = dd_CreateMatrix(0, dimension + 1);
  dd_SetMatrixRepresentationType(intersectionH, dd_Inequality);

  // for each convex
  for (UnsignedInteger i = 0; i < remainingSize; ++ i)
  {
    const Sample vertices1(coll[remainingIndices[i]]);
    const UnsignedInteger nv1 = vertices1.getSize();

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
    for (UnsignedInteger i = 0; i < intersectionVerticesNumber; ++i)
    {
      // First entry = 1 -> point, 0 -> ray
      if (dd_get_d(gen->matrix[i][0]) != 1.0)
        throw InternalException(HERE) << "assumed only points, no rays";

      Point vertex(dimension);
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        vertex[j] = dd_get_d(gen->matrix[i][j + 1]);
      result.add(vertex);
    }
  }
  dd_FreeMatrix(gen);

  return result;
#else
  throw NotYetImplementedException(HERE) << "No cddlib support";
#endif
}


// we intersect the convex cylinders in one pass first (if any) to initialize a list of unions of convexes represented as their vertices
// the main loop iterates the list of non-convex cylinders, with each cylinder is decomposed as a union of convexes
// if there were no convex cylinders at the previous step we decompose the first non-convex cylinder to initialize the list of unions
// then for each new cylinder we compute all the intersections of each convex in the current union list with each convex in its decomposition
// then the current list of unions is updated, removing the empty intersections
// finally the last union of convexes obtained after visiting all cylinders is assembled in a single mesh
Mesh IntersectionMesher::buildCylinder(const Collection<Cylinder> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));

  // intersect all convex cylinders first
  Collection<Sample> unionCurrent;
  Indices nonConvex;
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    if (coll[i].isConvex())
      unionCurrent.add(coll[i].getVertices());
    else
      nonConvex.add(i);
  }
  if (unionCurrent.getSize() > 1)
  {
    const Sample convexIntersection(buildConvexSample(unionCurrent));
    if (!convexIntersection.getSize())
      return Mesh(convexIntersection);
    unionCurrent = {convexIntersection};
  }

  // if there are only non-convex cylinders we must decompose the first one to initialize unionCurrent
  ConvexDecompositionMesher convexDecompositionMesher;
  convexDecompositionMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
  const UnsignedInteger nonConvexSize = nonConvex.getSize();
  UnsignedInteger startNonConvex = 0;
  if (nonConvexSize == size)
  {
    // build decomposition of first non-convex cylinder
    const Cylinder cylinder0(coll[nonConvex[0]]);
    const Collection<Mesh> baseDecomposition(convexDecompositionMesher.build(cylinder0.getBase()));
    const UnsignedInteger baseDecompositionSize = baseDecomposition.getSize();
    for (UnsignedInteger j = 0; j < baseDecompositionSize; ++ j)
    {
      const Cylinder cylinder0J(baseDecomposition[j],
                                cylinder0.getExtension(),
                                cylinder0.getInjection(),
                                cylinder0.getDiscretization());
      unionCurrent.add(cylinder0J.getVertices());
    }

    // start at index 1
    startNonConvex = 1;
  }

  // for each remaining (non-convex) cylinder i
  for (UnsignedInteger i = startNonConvex; i < nonConvexSize; ++ i)
  {
    // build decomposition
    Collection<Sample> unionNext;
    const Cylinder cylinderI(coll[nonConvex[i]]);
    const Collection<Mesh> baseDecomposition(convexDecompositionMesher.build(cylinderI.getBase()));
    const UnsignedInteger baseDecompositionSize = baseDecomposition.getSize();
    for (UnsignedInteger j = 0; j < baseDecompositionSize; ++ j)
    {
      const Cylinder cylinderIJ(baseDecomposition[j],
                                cylinderI.getExtension(),
                                cylinderI.getInjection(),
                                cylinderI.getDiscretization());
      unionNext.add(cylinderIJ.getVertices());
    }

    // build lists of intersections to compute
    const UnsignedInteger toDoSize = unionCurrent.getSize() *  unionNext.getSize();

    // loop over intersections
    Collection<Sample> result(toDoSize);
    const IntersectionMesherConvexSamplePolicy policy(*this, unionCurrent, unionNext, 0, result);
    TBBImplementation::ParallelFor(0, toDoSize, policy);

    // prune empty intersections
    unionCurrent.resize(0);
    for (UnsignedInteger i0 = 0; i0 < result.getSize(); ++ i0)
      if (result[i0].getSize())
        unionCurrent.add(result[i0]);

    // early exit if there are no non-empty intersections at this stage
    if (!unionCurrent.getSize())
      return Mesh(Sample(0, cylinderI.getDimension()));

  } // for cylinder i

  // build mesh of union of remaining intersections
  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(unionCurrent.getSize());
  for (UnsignedInteger i = 0; i < unionCurrent.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(unionCurrent[i]);
  return UnionMesher().build(collMesh);
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

/* Simplices decomposition flag accessor */
void IntersectionMesher::setUseSimplicesDecomposition(const Bool useSimplicesDecomposition)
{
  useSimplicesDecomposition_ = useSimplicesDecomposition;
}

Bool IntersectionMesher::getUseSimplicesDecomposition() const
{
  return useSimplicesDecomposition_;
}

/* Method save() stores the object through the StorageManager */
void IntersectionMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("recompress_", recompress_);
  adv.saveAttribute("useSimplicesDecomposition_", useSimplicesDecomposition_);
}

/* Method load() reloads the object from the StorageManager */
void IntersectionMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute("recompress_", recompress_);
  adv.loadAttribute("useSimplicesDecomposition_", useSimplicesDecomposition_);
}

struct IntersectionMesher_init
{
  IntersectionMesher_init()
  {
    ResourceMap::AddAsUnsignedInteger("IntersectionMesher-BlockSize", 1 << 16);
#ifdef OPENTURNS_HAVE_CDDLIB
    dd_set_global_constants();
#endif
  }

  ~IntersectionMesher_init()
  {
#ifdef OPENTURNS_HAVE_CDDLIB
    dd_free_global_constants();
#endif
  }
};

static IntersectionMesher_init __IntersectionMesher_initializer;

}
