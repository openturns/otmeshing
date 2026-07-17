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
#include <openturns/TBBImplementation.hxx>

#include "otmeshing/IntersectionMesher.hxx"
#include "otmeshing/CloudMesher.hxx"
#include "otmeshing/ConvexDecompositionMesher.hxx"
#include "otmeshing/ConvexHullMesher.hxx"
#include "otmeshing/UnionMesher.hxx"
#include "otmeshing/VolumeMesher.hxx"

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
      output_[n] = intersectionMesher_.buildConvexSample(input1_[i], input2_[j]);
    }
  }
};

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

struct ScopedMatrixPtr
{
  dd_MatrixPtr ptr_;

  ScopedMatrixPtr() : ptr_(nullptr) {}
  explicit ScopedMatrixPtr(dd_MatrixPtr p) : ptr_(p) {}
  ~ScopedMatrixPtr() { if (ptr_) dd_FreeMatrix(ptr_); }

  ScopedMatrixPtr(ScopedMatrixPtr && other) : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  ScopedMatrixPtr(const ScopedMatrixPtr &) = delete;
  ScopedMatrixPtr & operator=(ScopedMatrixPtr && other)
  {
    if (ptr_) dd_FreeMatrix(ptr_);
    ptr_ = other.ptr_;
    other.ptr_ = nullptr;
    return *this;
  }
  ScopedMatrixPtr & operator=(const ScopedMatrixPtr &) = delete;

  dd_MatrixPtr get() const { return ptr_; }
  dd_MatrixPtr * ptrAddr() { return &ptr_; }
  dd_MatrixPtr release() { dd_MatrixPtr p = ptr_; ptr_ = nullptr; return p; }
  void reset(dd_MatrixPtr p = nullptr) { if (ptr_) dd_FreeMatrix(ptr_); ptr_ = p; }
  dd_MatrixPtr operator->() const { return ptr_; }
};

struct ScopedPolyhedraPtr
{
  dd_PolyhedraPtr ptr_;

  ScopedPolyhedraPtr() : ptr_(nullptr) {}
  explicit ScopedPolyhedraPtr(dd_PolyhedraPtr p) : ptr_(p) {}
  ~ScopedPolyhedraPtr() { if (ptr_) dd_FreePolyhedra(ptr_); }

  ScopedPolyhedraPtr(ScopedPolyhedraPtr && other) : ptr_(other.ptr_) { other.ptr_ = nullptr; }
  ScopedPolyhedraPtr(const ScopedPolyhedraPtr &) = delete;
  ScopedPolyhedraPtr & operator=(ScopedPolyhedraPtr && other)
  {
    if (ptr_) dd_FreePolyhedra(ptr_);
    ptr_ = other.ptr_;
    other.ptr_ = nullptr;
    return *this;
  }
  ScopedPolyhedraPtr & operator=(const ScopedPolyhedraPtr &) = delete;

  dd_PolyhedraPtr get() const { return ptr_; }
  dd_PolyhedraPtr release() { dd_PolyhedraPtr p = ptr_; ptr_ = nullptr; return p; }
  void reset(dd_PolyhedraPtr p = nullptr) { if (ptr_) dd_FreePolyhedra(ptr_); ptr_ = p; }
  dd_PolyhedraPtr operator->() const { return ptr_; }
};

static dd_MatrixPtr ComputeHRepresentation(const Sample & vertices)
{
  const UnsignedInteger nv = vertices.getSize();
  const UnsignedInteger dimension = vertices.getDimension();
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr m(dd_CreateMatrix(nv, dimension + 1));
  dd_SetMatrixRepresentationType(m.get(), dd_Generator);
  for (UnsignedInteger i = 0; i < nv; ++ i)
  {
    dd_set_d(m->matrix[i][0], 1.0);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      dd_set_d(m->matrix[i][k + 1], vertices(i, k));
  }

  ScopedPolyhedraPtr p(dd_DDMatrix2Poly(m.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed: " << cdd_error_to_string(err);

  return dd_CopyInequalities(p.get());
}

static dd_MatrixPtr ComputeHRepresentation(const Sample & allVertices,
    const IndicesCollection & simplices,
    const UnsignedInteger simplexIndex)
{
  const UnsignedInteger dimension = allVertices.getDimension();
  const UnsignedInteger nv = dimension + 1;
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr m(dd_CreateMatrix(nv, dimension + 1));
  dd_SetMatrixRepresentationType(m.get(), dd_Generator);
  for (UnsignedInteger i = 0; i < nv; ++ i)
  {
    dd_set_d(m->matrix[i][0], 1.0);
    const UnsignedInteger vertexIndex = simplices(simplexIndex, i);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      dd_set_d(m->matrix[i][k + 1], allVertices(vertexIndex, k));
  }

  ScopedPolyhedraPtr p(dd_DDMatrix2Poly(m.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed: " << cdd_error_to_string(err);

  return dd_CopyInequalities(p.get());
}

static Sample IntersectFromH(dd_MatrixPtr h1, dd_MatrixPtr h2, const UnsignedInteger dimension)
{
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr intersectionH(dd_CreateMatrix(0, dimension + 1));
  dd_SetMatrixRepresentationType(intersectionH.get(), dd_Inequality);

  dd_MatrixAppendTo(intersectionH.ptrAddr(), h1);
  dd_MatrixAppendTo(intersectionH.ptrAddr(), h2);

  ScopedPolyhedraPtr intersectionV(dd_DDMatrix2Poly(intersectionH.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for intersection: " << cdd_error_to_string(err);

  ScopedMatrixPtr gen(dd_CopyGenerators(intersectionV.get()));

  Sample result(0, dimension);
  const UnsignedInteger intersectionVerticesNumber = gen->rowsize;
  if (intersectionVerticesNumber >= (dimension + 1))
  {
    result = Sample(intersectionVerticesNumber, dimension);
    for (UnsignedInteger i = 0; i < intersectionVerticesNumber; ++ i)
    {
      if (dd_get_d(gen->matrix[i][0]) != 1.0)
        throw InternalException(HERE) << "assumed only points, no rays";
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        result(i, j) = dd_get_d(gen->matrix[i][j + 1]);
    }
  }
  return result;
}

struct IntersectionMesherConvexSampleHMatrixPolicy
{
  const std::vector<ScopedMatrixPtr> & h1_;
  const std::vector<ScopedMatrixPtr> & h2_;
  UnsignedInteger done_;
  Collection<Sample> & output_;
  UnsignedInteger stride_;
  UnsignedInteger dimension_;

  IntersectionMesherConvexSampleHMatrixPolicy(const std::vector<ScopedMatrixPtr> & h1,
      const std::vector<ScopedMatrixPtr> & h2,
      const UnsignedInteger done,
      Collection<Sample> & output,
      const UnsignedInteger dimension)
    : h1_(h1)
    , h2_(h2)
    , done_(done)
    , output_(output)
    , stride_(h1.size())
    , dimension_(dimension)
  {}

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger n = r.begin(); n != r.end(); ++n)
    {
      const UnsignedInteger i = (n + done_) % stride_;
      const UnsignedInteger j = (n + done_) / stride_;
      output_[n] = IntersectFromH(h1_[i].get(), h2_[j].get(), dimension_);
    }
  }
};
#endif

Mesh IntersectionMesher::build(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

  const UnsignedInteger dimension = coll[0].getDimension();

  Collection<Sample> unionVertices;
  UnsignedInteger unionSize = 0;
#ifdef OPENTURNS_HAVE_CDDLIB
  std::vector<ScopedMatrixPtr> unionH;
#endif
  {
    const IndicesCollection & simplices0 = coll[0].getSimplices();
    unionSize = simplices0.getSize();
    const Sample & meshVerts0 = coll[0].getVertices();
#ifdef OPENTURNS_HAVE_CDDLIB
    unionH.reserve(unionSize);
    for (UnsignedInteger j = 0; j < unionSize; ++ j)
      unionH.emplace_back(ComputeHRepresentation(meshVerts0, simplices0, j));
#else
    const UnsignedInteger stride0 = dimension + 1;
    for (UnsignedInteger j = 0; j < unionSize; ++ j)
    {
      Sample verts(stride0, dimension);
      for (UnsignedInteger k = 0; k < stride0; ++ k)
      {
        const UnsignedInteger vi = simplices0(j, k);
        for (UnsignedInteger d = 0; d < dimension; ++ d)
          verts(k, d) = meshVerts0(vi, d);
      }
      unionVertices.add(verts);
    }
#endif
  }

  for (UnsignedInteger i = 1; i < size; ++i)
  {
    UnsignedInteger nextSize = 0;
#ifdef OPENTURNS_HAVE_CDDLIB
    std::vector<ScopedMatrixPtr> nextH;
#else
    Collection<Sample> nextVertices;
#endif
    {
      const IndicesCollection & simplicesI = coll[i].getSimplices();
      nextSize = simplicesI.getSize();
      const Sample & meshVertsI = coll[i].getVertices();
#ifdef OPENTURNS_HAVE_CDDLIB
      nextH.reserve(nextSize);
      for (UnsignedInteger j = 0; j < nextSize; ++ j)
        nextH.emplace_back(ComputeHRepresentation(meshVertsI, simplicesI, j));
#else
      const UnsignedInteger strideI = dimension + 1;
      for (UnsignedInteger j = 0; j < nextSize; ++ j)
      {
        Sample verts(strideI, dimension);
        for (UnsignedInteger k = 0; k < strideI; ++ k)
        {
          const UnsignedInteger vi = simplicesI(j, k);
          for (UnsignedInteger d = 0; d < dimension; ++ d)
            verts(k, d) = meshVertsI(vi, d);
        }
        nextVertices.add(verts);
      }
#endif
    }

    const UnsignedInteger blockSize = std::max(UnsignedInteger(1), ResourceMap::GetAsUnsignedInteger("IntersectionMesher-BlockSize"));
    Collection<Sample> result(0);
    const UnsignedInteger toDoSize = unionSize * nextSize;

    Collection<Sample> resultChunk(std::min(blockSize, toDoSize));
    for (UnsignedInteger done = 0; done < toDoSize; done += blockSize)
    {
      const UnsignedInteger actualBlockSize = std::min(blockSize, toDoSize - done);
      resultChunk.resize(actualBlockSize);
#ifdef OPENTURNS_HAVE_CDDLIB
      const IntersectionMesherConvexSampleHMatrixPolicy policy(unionH, nextH, done, resultChunk, dimension);
#else
      const IntersectionMesherConvexSamplePolicy policy(*this, unionVertices, nextVertices, done, resultChunk);
#endif
      TBBImplementation::ParallelFor(0, actualBlockSize, policy);

      for (UnsignedInteger i0 = 0; i0 < actualBlockSize; ++ i0)
      {
        Sample& sample = resultChunk[i0];
        if (sample.getSize())
          result.add(sample);
      }
    }

    unionVertices = result;
    unionSize = unionVertices.getSize();
    if (!unionSize)
      return Mesh(Sample(0, dimension));

#ifdef OPENTURNS_HAVE_CDDLIB
    unionH.clear();
    unionH.reserve(unionSize);
    for (UnsignedInteger k = 0; k < unionSize; ++ k)
      unionH.emplace_back(ComputeHRepresentation(unionVertices[k]));
#endif
  }

  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(unionVertices.getSize());
  for (UnsignedInteger i = 0; i < unionVertices.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(unionVertices[i]);
  return UnionMesher().build(collMesh);
}

Mesh IntersectionMesher::buildWithConvexParts(const Mesh & mesh, const SampleCollection & convexPieces) const
{
  ConvexDecompositionMesher convexDecompositionMesher;
  convexDecompositionMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
  Collection<Sample> unionCurrent;
  const Collection<Mesh> baseDecomposition0(convexDecompositionMesher.build(mesh));
  for (UnsignedInteger j = 0; j < baseDecomposition0.getSize(); ++ j)
    unionCurrent.add(baseDecomposition0[j].getVertices());

  const UnsignedInteger toDoSize = unionCurrent.getSize() * convexPieces.getSize();
  const UnsignedInteger blockSize = ResourceMap::GetAsUnsignedInteger("IntersectionMesher-BlockSize");
  Collection<Sample> result(0);
  Collection<Sample> resultChunk(blockSize);
  for (UnsignedInteger done = 0; done < toDoSize; done += blockSize)
  {
    const UnsignedInteger actualBlockSize = std::min(blockSize, toDoSize - done);
    const IntersectionMesherConvexSamplePolicy policy(*this, unionCurrent, convexPieces, done, resultChunk);
    TBBImplementation::ParallelFor(0, actualBlockSize, policy);
    for (UnsignedInteger i0 = 0; i0 < actualBlockSize; ++ i0)
      if (resultChunk[i0].getSize())
        result.add(resultChunk[i0]);
  }

  if (!result.getSize())
    return Mesh(Sample(0, mesh.getDimension()));

  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(result.getSize());
  for (UnsignedInteger i = 0; i < result.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(result[i]);
  return UnionMesher().build(collMesh);
}

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
    // fewer tets than delaunay with CloudMesher
    result = VolumeMesher().build(ConvexHullMesher().build(intersectionVertices));
  }
  return result;
}


Sample IntersectionMesher::buildConvexSample(const Sample & s1, const Sample & s2) const
{
  const UnsignedInteger dimension = s1.getDimension();
  if (s2.getDimension() != dimension)
    throw InvalidArgumentException(HERE) << "IntersectionMesher expected vertices of same dimension";

  Sample result(0, dimension);
  const UnsignedInteger nv1 = s1.getSize();
  const UnsignedInteger nv2 = s2.getSize();

  Point min1(dimension, SpecFunc::Infinity);
  Point max1(dimension, -SpecFunc::Infinity);
  for (UnsignedInteger i = 0; i < nv1; ++ i)
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      const Scalar val = s1(i, k);
      if (val < min1[k]) min1[k] = val;
      if (val > max1[k]) max1[k] = val;
    }

  Point min2(dimension, SpecFunc::Infinity);
  Point max2(dimension, -SpecFunc::Infinity);
  for (UnsignedInteger i = 0; i < nv2; ++ i)
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      const Scalar val = s2(i, k);
      if (val < min2[k]) min2[k] = val;
      if (val > max2[k]) max2[k] = val;
    }

  for (UnsignedInteger k = 0; k < dimension; ++ k)
    if (std::max(min1[k], min2[k]) >= std::min(max1[k], max2[k]))
      return result;

#ifdef OPENTURNS_HAVE_CDDLIB
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr m1(dd_CreateMatrix(nv1, dimension + 1));
  dd_SetMatrixRepresentationType(m1.get(), dd_Generator);
  for (UnsignedInteger i1 = 0; i1 < nv1; ++ i1)
  {
    dd_set_d(m1->matrix[i1][0], 1.0);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      dd_set_d(m1->matrix[i1][k + 1], s1(i1, k));
  }
  ScopedPolyhedraPtr p1(dd_DDMatrix2Poly(m1.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for convex 1: " << cdd_error_to_string(err);
  ScopedMatrixPtr h1(dd_CopyInequalities(p1.get()));

  ScopedMatrixPtr m2(dd_CreateMatrix(nv2, dimension + 1));
  dd_SetMatrixRepresentationType(m2.get(), dd_Generator);
  for (UnsignedInteger i2 = 0; i2 < nv2; ++ i2)
  {
    dd_set_d(m2->matrix[i2][0], 1.0);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      dd_set_d(m2->matrix[i2][k + 1], s2(i2, k));
  }
  ScopedPolyhedraPtr p2(dd_DDMatrix2Poly(m2.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for convex 2: " << cdd_error_to_string(err);
  ScopedMatrixPtr h2(dd_CopyInequalities(p2.get()));

  return IntersectFromH(h1.get(), h2.get(), dimension);
#else
  throw NotYetImplementedException(HERE) << "No cddlib support";
#endif
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

  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr intersectionH(dd_CreateMatrix(0, dimension + 1));
  dd_SetMatrixRepresentationType(intersectionH.get(), dd_Inequality);

  for (UnsignedInteger i = 0; i < remainingSize; ++ i)
  {
    const Sample vertices1(coll[remainingIndices[i]]);
    const UnsignedInteger nv1 = vertices1.getSize();

    ScopedMatrixPtr m1(dd_CreateMatrix(nv1, dimension + 1));
    dd_SetMatrixRepresentationType(m1.get(), dd_Generator);
    for (UnsignedInteger i1 = 0; i1 < nv1; ++ i1)
    {
      dd_set_d(m1->matrix[i1][0], 1.0);
      for (UnsignedInteger k = 0; k < dimension; ++ k)
        dd_set_d(m1->matrix[i1][k + 1], vertices1(i1, k));
    }

    ScopedPolyhedraPtr p1(dd_DDMatrix2Poly(m1.get(), &err));
    if (err != dd_NoError)
      throw InternalException(HERE) << "dd_DDMatrix2Poly failed for mesh 1: " << cdd_error_to_string(err);

    ScopedMatrixPtr h1(dd_CopyInequalities(p1.get()));
    dd_MatrixAppendTo(intersectionH.ptrAddr(), h1.get());
  }

  ScopedPolyhedraPtr intersectionV(dd_DDMatrix2Poly(intersectionH.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for intersection: " << cdd_error_to_string(err);

  ScopedMatrixPtr gen(dd_CopyGenerators(intersectionV.get()));
  const UnsignedInteger intersectionVerticesNumber = gen->rowsize;
  if (intersectionVerticesNumber >= (dimension + 1))
  {
    result = Sample(intersectionVerticesNumber, dimension);
    for (UnsignedInteger i = 0; i < intersectionVerticesNumber; ++ i)
    {
      if (dd_get_d(gen->matrix[i][0]) != 1.0)
        throw InternalException(HERE) << "assumed only points, no rays";
      for (UnsignedInteger j = 0; j < dimension; ++ j)
        result(i, j) = dd_get_d(gen->matrix[i][j + 1]);
    }
  }
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
// the last union of convexes obtained after visiting all cylinders is returned

// Core logic for cylinder intersection, returning the convex decomposition
Collection<Sample> IntersectionMesher::buildCylinderConvex(const Collection<Cylinder> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Collection<Sample>(0);

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
      return Collection<Sample>(0);
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

    for (UnsignedInteger j = 0; j < baseDecomposition.getSize(); ++ j)
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
    for (UnsignedInteger j = 0; j < baseDecomposition.getSize(); ++ j)
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
      return Collection<Sample>(0);
  } // for cylinder i

  return unionCurrent;
}

// intersect cylinders and assemble into a single mesh
Mesh IntersectionMesher::buildCylinder(const Collection<Cylinder> & coll) const
{
  const SampleCollection convexPieces = buildCylinderConvex(coll);
  if (!convexPieces.getSize())
    return Mesh(Sample(0, coll.getSize() ? coll[0].getDimension() : 0));

  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(convexPieces.getSize());
  for (UnsignedInteger i = 0; i < convexPieces.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(convexPieces[i]);
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
    ResourceMap::AddAsScalar("ConvexDecompositionMesher-Threshold", 0.05);
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
