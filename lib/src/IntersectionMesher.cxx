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
#include <unordered_map>
#include <unordered_set>
#include <set>

#include <openturns/IntervalMesher.hxx>
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

// RAII wrappers for cddlib C resources.
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
  dd_SetMatrixNumberType(m.get(), dd_Rational);
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

  dd_MatrixPtr h = dd_CopyInequalities(p.get());
  return h;
}

// Compute H-representation for a simplex defined by an index into an IndicesCollection.
// Avoids creating an intermediate Sample for the simplex vertices.
static dd_MatrixPtr ComputeHRepresentation(const Sample & allVertices, const IndicesCollection & simplices, const UnsignedInteger simplexIndex)
{
  const UnsignedInteger dimension = allVertices.getDimension();
  const UnsignedInteger nv = dimension + 1;
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr m(dd_CreateMatrix(nv, dimension + 1));
  dd_SetMatrixNumberType(m.get(), dd_Rational);
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

  dd_MatrixPtr h = dd_CopyInequalities(p.get());
  return h;
}

static Sample IntersectFromH(dd_MatrixPtr h1, dd_MatrixPtr h2, const UnsignedInteger dimension)
{
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr intersectionH(dd_CreateMatrix(0, dimension + 1));
  dd_SetMatrixNumberType(intersectionH.get(), dd_Rational);
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

// Filtered versions that iterate over a pre-computed list of candidate pairs
// instead of the full cartesian product.
#ifdef OPENTURNS_HAVE_CDDLIB
struct IntersectionMesherConvexSampleHMatrixFilteredPolicy
{
  const std::vector<ScopedMatrixPtr> & h1_;
  const std::vector<ScopedMatrixPtr> & h2_;
  const std::vector<UnsignedInteger> & candI_;
  const std::vector<UnsignedInteger> & candJ_;
  UnsignedInteger done_;
  Collection<Sample> & output_;
  UnsignedInteger dimension_;

  IntersectionMesherConvexSampleHMatrixFilteredPolicy(
      const std::vector<ScopedMatrixPtr> & h1,
      const std::vector<ScopedMatrixPtr> & h2,
      const std::vector<UnsignedInteger> & candI,
      const std::vector<UnsignedInteger> & candJ,
      const UnsignedInteger done,
      Collection<Sample> & output,
      const UnsignedInteger dimension)
    : h1_(h1)
    , h2_(h2)
    , candI_(candI)
    , candJ_(candJ)
    , done_(done)
    , output_(output)
    , dimension_(dimension)
  {}

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger n = r.begin(); n != r.end(); ++ n)
    {
      const UnsignedInteger i = candI_[n + done_];
      const UnsignedInteger j = candJ_[n + done_];
      output_[n] = IntersectFromH(h1_[i].get(), h2_[j].get(), dimension_);
    }
  }
};
#endif

struct IntersectionMesherConvexSampleFilteredPolicy
{
  const IntersectionMesher & intersectionMesher_;
  const Collection<Sample> & input1_;
  const Collection<Sample> & input2_;
  const std::vector<UnsignedInteger> & candI_;
  const std::vector<UnsignedInteger> & candJ_;
  UnsignedInteger done_;
  Collection<Sample> & output_;

  IntersectionMesherConvexSampleFilteredPolicy(const IntersectionMesher & intersectionMesher,
      const Collection<Sample> & input1,
      const Collection<Sample> & input2,
      const std::vector<UnsignedInteger> & candI,
      const std::vector<UnsignedInteger> & candJ,
      const UnsignedInteger done,
      Collection<Sample> & output)
    : intersectionMesher_(intersectionMesher)
    , input1_(input1)
    , input2_(input2)
    , candI_(candI)
    , candJ_(candJ)
    , done_(done)
    , output_(output)
  {}

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger n = r.begin(); n != r.end(); ++ n)
    {
      const UnsignedInteger i = candI_[n + done_];
      const UnsignedInteger j = candJ_[n + done_];
      output_[n] = intersectionMesher_.buildConvexSample(input1_[i], input2_[j]);
    }
  }
};

static UnsignedInteger ComputeFactorial(const UnsignedInteger n)
{
  UnsignedInteger result = 1;
  for (UnsignedInteger i = 2; i <= n; ++ i) result *= i;
  return result;
}

// -----------------------------------------------------------------------
// AABB spatial index for pruning simplex-simplex intersection pairs.
// -----------------------------------------------------------------------

// Axis-aligned bounding box for a simplex.
struct SimplexAABB
{
  OT::Point lower;
  OT::Point upper;
};

// Compute the AABB for a single simplex by reading vertex coordinates
// directly from the mesh vertex array via the simplex index.
static SimplexAABB ComputeSimplexAABB(const Sample & vertices,
    const IndicesCollection & simplices,
    const UnsignedInteger simplexIndex)
{
  const UnsignedInteger dim = vertices.getDimension();
  SimplexAABB aabb;
  aabb.lower = Point(dim, SpecFunc::Infinity);
  aabb.upper = Point(dim, -SpecFunc::Infinity);
  for (UnsignedInteger i = 0; i <= dim; ++ i)
  {
    const UnsignedInteger vi = simplices(simplexIndex, i);
    for (UnsignedInteger k = 0; k < dim; ++ k)
    {
      const Scalar val = vertices(vi, k);
      if (val < aabb.lower[k]) aabb.lower[k] = val;
      if (val > aabb.upper[k]) aabb.upper[k] = val;
    }
  }
  return aabb;
}

// Check whether two AABBs overlap (non-empty interior intersection).
static Bool AABBOverlap(const SimplexAABB & a, const SimplexAABB & b)
{
  const UnsignedInteger dim = a.lower.getDimension();
  for (UnsignedInteger k = 0; k < dim; ++ k)
    if (!(a.upper[k] > b.lower[k]) || !(b.upper[k] > a.lower[k]))
      return false;
  return true;
}

// Uniform spatial grid index over a set of SimplexAABBs.
// Used to quickly find candidate simplices that may intersect a query AABB.
class SimplexGrid
{
public:
  SimplexGrid(const std::vector<SimplexAABB> & aabbs, UnsignedInteger targetCellsPerDim)
    : dimension_(aabbs.empty() ? 1 : aabbs[0].lower.getDimension())
  {
    if (aabbs.empty()) return;

    // Overall bounding box of all AABBs
    Point globalMin(dimension_, SpecFunc::Infinity);
    Point globalMax(dimension_, -SpecFunc::Infinity);
    for (const auto & aabb : aabbs)
    {
      for (UnsignedInteger k = 0; k < dimension_; ++ k)
      {
        if (aabb.lower[k] < globalMin[k]) globalMin[k] = aabb.lower[k];
        if (aabb.upper[k] > globalMax[k]) globalMax[k] = aabb.upper[k];
      }
    }

    origin_ = globalMin;
    Scalar maxExtent = 0;
    for (UnsignedInteger k = 0; k < dimension_; ++ k)
    {
      const Scalar ext = globalMax[k] - globalMin[k];
      if (ext > maxExtent) maxExtent = ext;
    }
    // Guard against degenerate domain
    if (maxExtent <= 0.0) maxExtent = 1.0;

    cellSize_ = maxExtent / static_cast<Scalar>(targetCellsPerDim);
    if (cellSize_ <= 0.0) cellSize_ = 1.0;

    nCells_.resize(dimension_);
    strides_.resize(dimension_);
    UnsignedInteger stride = 1;
    for (UnsignedInteger k = 0; k < dimension_; ++ k)
    {
      nCells_[k] = static_cast<UnsignedInteger>(
          std::max(1.0, ceil((globalMax[k] - globalMin[k]) / cellSize_)));
      strides_[k] = stride;
      stride *= nCells_[k];
    }
    totalCells_ = stride;

    // Insert each AABB into all cells it overlaps
    for (UnsignedInteger idx = 0; idx < aabbs.size(); ++ idx)
    {
      const auto & aabb = aabbs[idx];
      Indices lo(dimension_), hi(dimension_);
      for (UnsignedInteger k = 0; k < dimension_; ++ k)
      {
        lo[k] = static_cast<UnsignedInteger>(std::max(0.0,
            floor((aabb.lower[k] - origin_[k]) / cellSize_)));
        hi[k] = static_cast<UnsignedInteger>(std::min(
            static_cast<double>(nCells_[k] - 1),
            floor((aabb.upper[k] - origin_[k]) / cellSize_)));
        if (hi[k] < lo[k]) hi[k] = lo[k];
      }

      // Visit all cells in the box [lo, hi]
      UnsignedInteger nBox = 1;
      for (UnsignedInteger k = 0; k < dimension_; ++ k)
        nBox *= (hi[k] - lo[k] + 1);

      for (UnsignedInteger c = 0; c < nBox; ++ c)
      {
        Indices coord(dimension_);
        UnsignedInteger tmp = c;
        UnsignedInteger linIdx = 0;
        for (UnsignedInteger k = 0; k < dimension_; ++ k)
        {
          const UnsignedInteger span = hi[k] - lo[k] + 1;
          coord[k] = lo[k] + (tmp % span);
          tmp /= span;
          linIdx += coord[k] * strides_[k];
        }
        cells_[linIdx].push_back(idx);
      }
    }
  }

  // Query: for each simplex index whose cell overlaps the query AABB,
  // call the provided functor.
  template <typename Func>
  void query(const SimplexAABB & aabb, Func func) const
  {
    Indices lo(dimension_), hi(dimension_);
    for (UnsignedInteger k = 0; k < dimension_; ++ k)
    {
      const SignedInteger lo_s = std::max(SignedInteger(0),
          (SignedInteger)std::floor((aabb.lower[k] - origin_[k]) / cellSize_));
      const SignedInteger hi_s = std::min(
          (SignedInteger)(nCells_[k] - 1),
          (SignedInteger)std::floor((aabb.upper[k] - origin_[k]) / cellSize_));
      if (hi_s < lo_s) return; // no overlapping cell — AABB outside grid
      lo[k] = (UnsignedInteger)lo_s;
      hi[k] = (UnsignedInteger)hi_s;
    }

    UnsignedInteger nBox = 1;
    for (UnsignedInteger k = 0; k < dimension_; ++ k)
      nBox *= (hi[k] - lo[k] + 1);

    for (UnsignedInteger c = 0; c < nBox; ++ c)
    {
      Indices coord(dimension_);
      UnsignedInteger tmp = c;
      UnsignedInteger linIdx = 0;
      for (UnsignedInteger k = 0; k < dimension_; ++ k)
      {
        const UnsignedInteger span = hi[k] - lo[k] + 1;
        coord[k] = lo[k] + (tmp % span);
        tmp /= span;
        linIdx += coord[k] * strides_[k];
      }
      auto it = cells_.find(linIdx);
      if (it != cells_.end())
      {
        for (UnsignedInteger j : it->second)
          func(j);
      }
    }
  }

private:
  UnsignedInteger dimension_;
  Point origin_;
  Scalar cellSize_;
  Indices nCells_;
  Indices strides_;
  UnsignedInteger totalCells_;
  std::unordered_map<UnsignedInteger, std::vector<UnsignedInteger>> cells_;
};

// Return true if the mesh is a tensor-product grid produced by IntervalMesher.
// Validates that the vertex layout follows the lexicographic ordering
// convention (dim 0 varying fastest) with constant stride, and that the
// simplex count matches the implied grid resolution.
static Bool IsTensorProductGrid(const Mesh & mesh)
{
  const UnsignedInteger dimension = mesh.getDimension();
  const Sample & verts = mesh.getVertices();
  const UnsignedInteger vertexCount = verts.getSize();
  const UnsignedInteger simplicesNumber = mesh.getSimplicesNumber();
  const UnsignedInteger factor = ComputeFactorial(dimension);

  if (simplicesNumber % factor != 0)
    return false;

  const Point lower = verts.getMin();
  const Point upper = verts.getMax();

  Indices nk(dimension);
  UnsignedInteger stride = 1;
  for (UnsignedInteger k = 0; k < dimension; ++ k)
  {
    if (stride >= vertexCount)
      return false;
    const Scalar stepK = verts(stride, k) - verts(0, k);
    if (!(stepK > 0.0))
      return false;
    nk[k] = std::max(static_cast<UnsignedInteger>(1),
        static_cast<UnsignedInteger>(round((upper[k] - lower[k]) / stepK)));
    stride *= (nk[k] + 1);
  }
  if (stride != vertexCount)
    return false;

  // Verify simplex count matches the inferred grid: dimension! * product(nk)
  UnsignedInteger expectedSimplices = factor;
  for (UnsignedInteger k = 0; k < dimension; ++ k)
    expectedSimplices *= nk[k];
  return simplicesNumber == expectedSimplices;
}

// Fast path for IntervalMesher-derived meshes (full hyper-rectangles).
// Their intersection is simply the bounding-box overlap — no per-simplex
// cdd intersection required.
static Bool TryBuildIntervalIntersection(const Collection<Mesh> & coll,
    Mesh & result)
{
  const UnsignedInteger size = coll.getSize();
  const UnsignedInteger dimension = coll[0].getDimension();

  // Validate that every input mesh is a genuine tensor-product grid.
  for (UnsignedInteger i = 0; i < size; ++ i)
    if (!IsTensorProductGrid(coll[i]))
      return false;

  // Compute the intersection bounding box across all meshes.
  // Since each mesh is a full hyper-rectangle, the intersection is simply
  // the element-wise max of the lower bounds and min of the upper bounds.
  const Sample & verts0 = coll[0].getVertices();
  Point interLower = verts0.getMin();
  Point interUpper = verts0.getMax();
  for (UnsignedInteger i = 1; i < size; ++ i)
  {
    const Sample & verts = coll[i].getVertices();
    const Point lower = verts.getMin();
    const Point upper = verts.getMax();
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      if (lower[k] > interLower[k]) interLower[k] = lower[k];
      if (upper[k] < interUpper[k]) interUpper[k] = upper[k];
    }
  }

  // Empty intersection?
  for (UnsignedInteger k = 0; k < dimension; ++ k)
    if (!(interLower[k] < interUpper[k]))
    {
      result = Mesh(Sample(0, dimension));
      return true;
    }

  // Compute original grid resolution from mesh 0 vertex layout.
  // IntervalMesher stores vertices lexicographically, dim 0 varying fastest.
  // stride_k = product of (nk_i + 1) for i < k
  // step_k = verts0(stride_k, k) - verts0(0, k)
  const UnsignedInteger vertexCount = verts0.getSize();
  const Point lower0 = verts0.getMin();
  const Point upper0 = verts0.getMax();
  Indices nkOriginal(dimension);
  UnsignedInteger stride = 1;
  for (UnsignedInteger k = 0; k < dimension; ++ k)
  {
    if (stride >= vertexCount)
      return false;
    const Scalar stepK = verts0(stride, k) - verts0(0, k);
    if (!(stepK > 0.0))
      return false;
    nkOriginal[k] = std::max(static_cast<UnsignedInteger>(1),
        static_cast<UnsignedInteger>(
            round((upper0[k] - lower0[k]) / stepK)));
    stride *= (nkOriginal[k] + 1);
  }
  if (stride != vertexCount)
    return false;

  // Compute result resolution from intersection width and original step
  Indices resolution(dimension);
  for (UnsignedInteger k = 0; k < dimension; ++ k)
  {
    const Scalar stepK = (upper0[k] - lower0[k]) / nkOriginal[k];
    resolution[k] = std::max(static_cast<UnsignedInteger>(1),
        static_cast<UnsignedInteger>(
            round((interUpper[k] - interLower[k]) / stepK)));
  }

  // Build result via IntervalMesher
  IntervalMesher intervalMesher(resolution);
  result = intervalMesher.build(Interval(interLower, interUpper));
  return true;
}

Mesh IntersectionMesher::build(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

  const UnsignedInteger dimension = coll[0].getDimension();

  // Fast path for IntervalMesher-derived meshes
  {
    Mesh result;
    if (TryBuildIntervalIntersection(coll, result))
      return result;
  }

  // Build pieces from mesh 0.
  // When useSimplicesDecomposition_ is true we iterate over the mesh's own
  // IndicesCollection directly, avoiding intermediate Mesh object creation.
  // Otherwise we fall back to ConvexDecompositionMesher.
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

  // for each remaining mesh i
  for (UnsignedInteger i = 1; i < size; ++i)
  {
    // build pieces from mesh i's simplices
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

    // Build AABB spatial index for intersection pairs.
    // When both sets are large, precompute AABBs and use a uniform grid
    // to prune non-overlapping simplex pairs, avoiding the O(n*m) dense loop.
    const UnsignedInteger blockSize = std::max(UnsignedInteger(1), ResourceMap::GetAsUnsignedInteger("IntersectionMesher-BlockSize"));
    Collection<Sample> result(0);
    std::vector<UnsignedInteger> candI, candJ;
    UnsignedInteger toDoSize = unionSize * nextSize;
    static const UnsignedInteger SPATIAL_THRESHOLD = 50;
    Bool filtered = false;

    if ((unionSize > SPATIAL_THRESHOLD) && (nextSize > SPATIAL_THRESHOLD))
    {
      filtered = true;
      // Precompute AABBs for union pieces
      std::vector<SimplexAABB> unionAABB(unionSize);
#ifdef OPENTURNS_HAVE_CDDLIB
      // First iteration: union comes from mesh 0 simplices
      if (i == 1)
      {
        const IndicesCollection & simp0 = coll[0].getSimplices();
        const Sample & verts0 = coll[0].getVertices();
        for (UnsignedInteger j = 0; j < unionSize; ++ j)
          unionAABB[j] = ComputeSimplexAABB(verts0, simp0, j);
      }
      else
#endif
      {
        for (UnsignedInteger j = 0; j < unionSize; ++ j)
        {
          unionAABB[j].lower = unionVertices[j].getMin();
          unionAABB[j].upper = unionVertices[j].getMax();
        }
      }
      // Precompute AABBs for next pieces
      std::vector<SimplexAABB> nextAABB(nextSize);
#ifdef OPENTURNS_HAVE_CDDLIB
      {
        const IndicesCollection & simpI = coll[i].getSimplices();
        const Sample & vertsI = coll[i].getVertices();
        for (UnsignedInteger j = 0; j < nextSize; ++ j)
          nextAABB[j] = ComputeSimplexAABB(vertsI, simpI, j);
      }
#else
      for (UnsignedInteger j = 0; j < nextSize; ++ j)
      {
        nextAABB[j].lower = nextVertices[j].getMin();
        nextAABB[j].upper = nextVertices[j].getMax();
      }
#endif
      // Determine grid resolution: ~ (nextSize)^(1/d) cells per dim
      UnsignedInteger targetPerDim = static_cast<UnsignedInteger>(
          ceil(pow(static_cast<double>(nextSize), 1.0 / dimension)));
      if (targetPerDim < 4) targetPerDim = 4;
      if (targetPerDim > 25) targetPerDim = 25;

      // Heuristic reservation to reduce reallocation
      candI.reserve(std::min(unionSize * nextSize, unionSize * 128u));
      candJ.reserve(std::min(unionSize * nextSize, unionSize * 128u));

      SimplexGrid grid(nextAABB, targetPerDim);

      // Use a set to deduplicate — same pair may be emitted from multiple grid cells
      std::set<std::pair<UnsignedInteger, UnsignedInteger>> seen;
      for (UnsignedInteger j = 0; j < unionSize; ++ j)
      {
        grid.query(unionAABB[j], [&](UnsignedInteger k) {
          if (AABBOverlap(unionAABB[j], nextAABB[k]) && seen.insert({j, k}).second)
          {
            candI.push_back(j);
            candJ.push_back(k);
          }
        });
      }
    }
    if (filtered)
      toDoSize = candI.size();

    // Process (candidate) pairs in blocks
    Collection<Sample> resultChunk(std::min(blockSize, toDoSize));
    for (UnsignedInteger done = 0; done < toDoSize; done += blockSize)
    {
      const UnsignedInteger actualBlockSize = std::min(blockSize, toDoSize - done);
      resultChunk.resize(actualBlockSize);
      if (!filtered)
      {
        // Dense (no spatial pre-filtering) — use original (i,j) from modular arithmetic
#ifdef OPENTURNS_HAVE_CDDLIB
        const IntersectionMesherConvexSampleHMatrixPolicy policy(unionH, nextH, done, resultChunk, dimension);
#else
        const IntersectionMesherConvexSamplePolicy policy(*this, unionVertices, nextVertices, done, resultChunk);
#endif
        TBBImplementation::ParallelFor(0, actualBlockSize, policy);
      }
      else
      {
        // Filtered candidate-list path
#ifdef OPENTURNS_HAVE_CDDLIB
        const IntersectionMesherConvexSampleHMatrixFilteredPolicy policy(unionH, nextH, candI, candJ, done, resultChunk, dimension);
#else
        const IntersectionMesherConvexSampleFilteredPolicy policy(*this, unionVertices, nextVertices, candI, candJ, done, resultChunk);
#endif
        TBBImplementation::ParallelFor(0, actualBlockSize, policy);
      }

      // prune empty intersections
      for (UnsignedInteger i0 = 0; i0 < actualBlockSize; ++ i0)
      {
        Sample& sample = resultChunk[i0];
        if (sample.getSize())
          result.add(sample);
      }
    }

    // old unionH/nextH freed automatically via ScopedMatrixPtr destructors

    // early exit if there are no non-empty intersections at this stage
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
  } // for mesh i

  // build mesh of union of remaining intersections
  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(unionVertices.getSize());
  for (UnsignedInteger i = 0; i < unionVertices.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(unionVertices[i]);
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
    // V>d+1, decompose into several simplices
    result = CloudMesher().build(intersectionVertices);
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

  // Compute bbox of s1 and s2 in a single pass each
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

  // bbox pruning
  for (UnsignedInteger k = 0; k < dimension; ++ k)
    if (std::max(min1[k], min2[k]) >= std::min(max1[k], max2[k]))
      return result;

#ifdef OPENTURNS_HAVE_CDDLIB
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr intersectionH(dd_CreateMatrix(0, dimension + 1));
  dd_SetMatrixNumberType(intersectionH.get(), dd_Rational);
  dd_SetMatrixRepresentationType(intersectionH.get(), dd_Inequality);

  // Process s1
  {
    ScopedMatrixPtr m1(dd_CreateMatrix(nv1, dimension + 1));
    dd_SetMatrixNumberType(m1.get(), dd_Rational);
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
    dd_MatrixAppendTo(intersectionH.ptrAddr(), h1.get());
  }

  // Process s2
  {
    ScopedMatrixPtr m2(dd_CreateMatrix(nv2, dimension + 1));
    dd_SetMatrixNumberType(m2.get(), dd_Rational);
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
    dd_MatrixAppendTo(intersectionH.ptrAddr(), h2.get());
  }

  // Convert intersection back to V-representation
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
    const Sample& vertices1 = coll[i];
    const UnsignedInteger nv1 = vertices1.getSize();
    Point localMin(dimension, SpecFunc::Infinity);
    Point localMax(dimension, -SpecFunc::Infinity);
    for (UnsignedInteger i1 = 0; i1 < nv1; ++ i1)
    {
      for (UnsignedInteger k = 0; k < dimension; ++ k)
      {
        const Scalar val = vertices1(i1, k);
        if (val < localMin[k]) localMin[k] = val;
        if (val > localMax[k]) localMax[k] = val;
      }
    }
    Bool pruned = false;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      pruned = std::max(lower1[k], localMin[k]) >= std::min(upper1[k], localMax[k]);
      if (pruned)
        break;
    }
    if (pruned)
      continue;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      lower1[k] = std::max(lower1[k], localMin[k]);
      upper1[k] = std::min(upper1[k], localMax[k]);
    }
    remainingIndices.add(i);
  }

  const UnsignedInteger remainingSize = remainingIndices.getSize();
  if (remainingSize < size)
    return result;

#ifdef OPENTURNS_HAVE_CDDLIB

  // initialize cddlib
  dd_ErrorType err = dd_NoError;

  ScopedMatrixPtr intersectionH(dd_CreateMatrix(0, dimension + 1));
  dd_SetMatrixNumberType(intersectionH.get(), dd_Rational);
  dd_SetMatrixRepresentationType(intersectionH.get(), dd_Inequality);

  // for each convex
  for (UnsignedInteger i = 0; i < remainingSize; ++ i)
  {
    const Sample& vertices1 = coll[remainingIndices[i]];
    const UnsignedInteger nv1 = vertices1.getSize();

    ScopedMatrixPtr m1(dd_CreateMatrix(nv1, dimension + 1));
    dd_SetMatrixNumberType(m1.get(), dd_Rational);
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
  } // i loop

  ScopedPolyhedraPtr intersectionV(dd_DDMatrix2Poly(intersectionH.get(), &err));
  if (err != dd_NoError)
    throw InternalException(HERE) << "dd_DDMatrix2Poly failed for intersection: " << cdd_error_to_string(err);

  ScopedMatrixPtr gen(dd_CopyGenerators(intersectionV.get()));
  const UnsignedInteger intersectionVerticesNumber = gen->rowsize;
  if (intersectionVerticesNumber >= (dimension + 1))
  {
    result = Sample(intersectionVerticesNumber, dimension);
    for (UnsignedInteger i = 0; i < intersectionVerticesNumber; ++i)
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
// finally the last union of convexes obtained after visiting all cylinders is assembled in a single mesh
Mesh IntersectionMesher::buildCylinder(const Collection<Cylinder> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));

  // Fast path: when all cylinders share the same injection pattern and are convex,
  // the intersection decomposes into independent base and interval intersections:
  //   ∩_i (Base_i × Ext) = (∩_i Base_i) × Ext
  // The base intersection is solved in the lower base dimension, then the result
  // is extruded back to the full dimension.
  if (size >= 2)
  {
    const Indices & inj0 = coll[0].getInjection();
    const Interval & ext0 = coll[0].getExtension();
    const UnsignedInteger disc0 = coll[0].getDiscretization();
    Bool allSame = coll[0].isConvex();
    for (UnsignedInteger i = 1; i < size; ++ i)
    {
      if (!coll[i].isConvex()
          || !(coll[i].getInjection() == inj0)
          || !(coll[i].getExtension() == ext0)
          || coll[i].getDiscretization() != disc0)
      {
        allSame = false;
        break;
      }
    }
    if (allSame)
    {
      // Intersect the bases in their native (lower) dimension
      Collection<Sample> baseSamples(size);
      for (UnsignedInteger i = 0; i < size; ++ i)
        baseSamples[i] = coll[i].getBase().getVertices();
      const Sample baseIntersection(buildConvexSample(baseSamples));
      if (baseIntersection.getSize() == 0)
        return Mesh(Sample(0, coll[0].getDimension()));

      // Extrude the intersected base and build the full-dimension mesh
      const Mesh baseMesh = CloudMesher().build(baseIntersection);
      const Cylinder resultCyl(baseMesh, ext0, inj0, disc0);
      return CloudMesher().build(resultCyl.getVertices());
    }
  }

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

  // pre-decompose all non-convex cylinders' bases in parallel
  const UnsignedInteger nonConvexSize = nonConvex.getSize();
  Collection<Collection<Mesh>> nonConvexDecompositions(nonConvexSize);
  {
    ConvexDecompositionMesher convexDecompositionMesher;
    convexDecompositionMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
    TBBImplementation::ParallelFor(0, nonConvexSize, [&](const TBBImplementation::BlockedRange<UnsignedInteger> & r)
    {
      ConvexDecompositionMesher localMesher;
      localMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
      for (UnsignedInteger i = r.begin(); i != r.end(); ++ i)
      {
        const Cylinder cylinderI(coll[nonConvex[i]]);
        nonConvexDecompositions[i] = localMesher.build(cylinderI.getBase());
      }
    });
  }

  // Initialize H-representation cache for unionCurrent
#ifdef OPENTURNS_HAVE_CDDLIB
  std::vector<ScopedMatrixPtr> unionCurrentH;
  {
    const UnsignedInteger sz = unionCurrent.getSize();
    unionCurrentH.reserve(sz);
    for (UnsignedInteger k = 0; k < sz; ++ k)
      unionCurrentH.emplace_back(ComputeHRepresentation(unionCurrent[k]));
  }
#endif

  // if there are only non-convex cylinders we must decompose the first one to initialize unionCurrent
  UnsignedInteger startNonConvex = 0;
  if (nonConvexSize == size)
  {
    // build decomposition of first non-convex cylinder from pre-computed data
    const Cylinder cylinder0(coll[nonConvex[0]]);
    const Collection<Mesh>& baseDecomposition0 = nonConvexDecompositions[0];
    const UnsignedInteger baseDecompositionSize0 = baseDecomposition0.getSize();
    for (UnsignedInteger j = 0; j < baseDecompositionSize0; ++ j)
    {
      const Cylinder cylinder0J(baseDecomposition0[j],
                                cylinder0.getExtension(),
                                cylinder0.getInjection(),
                                cylinder0.getDiscretization());
      unionCurrent.add(cylinder0J.getVertices());
#ifdef OPENTURNS_HAVE_CDDLIB
      unionCurrentH.emplace_back(ComputeHRepresentation(cylinder0J.getVertices()));
#endif
    }

    // start at index 1
    startNonConvex = 1;
  }

  // for each remaining (non-convex) cylinder i
  for (UnsignedInteger i = startNonConvex; i < nonConvexSize; ++ i)
  {
    // build decomposition from pre-computed data
    Collection<Sample> unionNext;
#ifdef OPENTURNS_HAVE_CDDLIB
    std::vector<ScopedMatrixPtr> nextH;
#endif
    const Cylinder cylinderI(coll[nonConvex[i]]);
    const Collection<Mesh>& baseDecomposition = nonConvexDecompositions[i];
    const UnsignedInteger baseDecompositionSize = baseDecomposition.getSize();
#ifdef OPENTURNS_HAVE_CDDLIB
    nextH.reserve(baseDecompositionSize);
#endif
    for (UnsignedInteger j = 0; j < baseDecompositionSize; ++ j)
    {
      const Cylinder cylinderIJ(baseDecomposition[j],
                                cylinderI.getExtension(),
                                cylinderI.getInjection(),
                                cylinderI.getDiscretization());
      const Sample & verts = cylinderIJ.getVertices();
      unionNext.add(verts);
#ifdef OPENTURNS_HAVE_CDDLIB
      nextH.emplace_back(ComputeHRepresentation(verts));
#endif
    }

    // Build AABB spatial index for intersection pairs (same approach as build())
    const UnsignedInteger blockSize = std::max(UnsignedInteger(1), ResourceMap::GetAsUnsignedInteger("IntersectionMesher-BlockSize"));
    const UnsignedInteger toDoSize = unionCurrent.getSize() * unionNext.getSize();
    const UnsignedInteger dimension = cylinderI.getDimension();
    Collection<Sample> result(0);
    std::vector<UnsignedInteger> candI, candJ;
    static const UnsignedInteger SPATIAL_THRESHOLD = 50;
    Bool filtered = false;

    if ((unionCurrent.getSize() > SPATIAL_THRESHOLD) && (unionNext.getSize() > SPATIAL_THRESHOLD))
    {
      filtered = true;
      // Precompute AABBs for current pieces
      // We use getMin()/getMax() on each vertex Sample because cylinder
      // decomposition gives us vertex sets, not mesh+simplices.
      std::vector<SimplexAABB> currentAABB(unionCurrent.getSize());
      for (UnsignedInteger j = 0; j < unionCurrent.getSize(); ++ j)
      {
        currentAABB[j].lower = unionCurrent[j].getMin();
        currentAABB[j].upper = unionCurrent[j].getMax();
      }

      // Precompute AABBs for next pieces
      std::vector<SimplexAABB> nextAABB(unionNext.getSize());
      for (UnsignedInteger j = 0; j < unionNext.getSize(); ++ j)
      {
        nextAABB[j].lower = unionNext[j].getMin();
        nextAABB[j].upper = unionNext[j].getMax();
      }

      // Determine grid resolution
      UnsignedInteger targetPerDim = static_cast<UnsignedInteger>(
          ceil(pow(static_cast<double>(unionNext.getSize()), 1.0 / dimension)));
      if (targetPerDim < 4) targetPerDim = 4;
      if (targetPerDim > 25) targetPerDim = 25;

      candI.reserve(std::min(toDoSize, unionCurrent.getSize() * 128u));
      candJ.reserve(std::min(toDoSize, unionCurrent.getSize() * 128u));

      SimplexGrid grid(nextAABB, targetPerDim);

      std::set<std::pair<UnsignedInteger, UnsignedInteger>> seen;
      for (UnsignedInteger j = 0; j < unionCurrent.getSize(); ++ j)
      {
        grid.query(currentAABB[j], [&](UnsignedInteger k) {
          if (AABBOverlap(currentAABB[j], nextAABB[k]) && seen.insert({j, k}).second)
          {
            candI.push_back(j);
            candJ.push_back(k);
          }
        });
      }
    }

    const UnsignedInteger effectiveToDo = filtered ? candI.size() : toDoSize;
    Collection<Sample> resultChunk(std::min(blockSize, effectiveToDo));
    for (UnsignedInteger done = 0; done < effectiveToDo; done += blockSize)
    {
      const UnsignedInteger actualBlockSize = std::min(blockSize, effectiveToDo - done);
      resultChunk.resize(actualBlockSize);
      if (!filtered)
      {
#ifdef OPENTURNS_HAVE_CDDLIB
        const IntersectionMesherConvexSampleHMatrixPolicy policy(unionCurrentH, nextH, done, resultChunk, dimension);
#else
        const IntersectionMesherConvexSamplePolicy policy(*this, unionCurrent, unionNext, done, resultChunk);
#endif
        TBBImplementation::ParallelFor(0, actualBlockSize, policy);
      }
      else
      {
#ifdef OPENTURNS_HAVE_CDDLIB
        const IntersectionMesherConvexSampleHMatrixFilteredPolicy policy(unionCurrentH, nextH, candI, candJ, done, resultChunk, dimension);
#else
        const IntersectionMesherConvexSampleFilteredPolicy policy(*this, unionCurrent, unionNext, candI, candJ, done, resultChunk);
#endif
        TBBImplementation::ParallelFor(0, actualBlockSize, policy);
      }

      // prune empty intersections
      for (UnsignedInteger i0 = 0; i0 < actualBlockSize; ++ i0)
      {
        Sample & sample = resultChunk[i0];
        if (sample.getSize())
          result.add(sample);
      }
    }

    // old unionCurrentH/nextH freed automatically via ScopedMatrixPtr destructors

    // early exit if there are no non-empty intersections at this stage
    unionCurrent = result;
    if (!unionCurrent.getSize())
      return Mesh(Sample(0, dimension));

    // Pre-compute H for the new unionCurrent pieces for next iteration
#ifdef OPENTURNS_HAVE_CDDLIB
    unionCurrentH.clear();
    {
      const UnsignedInteger sz = unionCurrent.getSize();
      unionCurrentH.reserve(sz);
      for (UnsignedInteger k = 0; k < sz; ++ k)
        unionCurrentH.emplace_back(ComputeHRepresentation(unionCurrent[k]));
    }
#endif

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
