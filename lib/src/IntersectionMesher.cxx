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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

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

// Collapse a convex vertex cloud to its extreme points. Grid-like inputs
// carry masses of interior vertices that are redundant for the H-rep but
// make the exact V-H conversion (GMP rationals) explode; the hull is the
// same polytope, so the intersection is unchanged.
static Sample ReduceConvexCloud(const Sample & cloud)
{
  const UnsignedInteger dimension = cloud.getDimension();
  if (cloud.getSize() <= 2 * (dimension + 1))
    return cloud;
  // Only full-dimensional clouds: flat cells hang Qhull instead of failing
  const Point cellMin(cloud.getMin());
  const Point cellMax(cloud.getMax());
  for (UnsignedInteger k = 0; k < dimension; ++ k)
    if (!(cellMax[k] > cellMin[k]))
      return cloud;
  // Explicit hulls pay off only where Qhull is cheap (measured: instant in
  // low dim, 11 s for 64 points in 10D, 329 s for 128; full 1024-point 10D
  // cloud hangs). Above dim 8 the exact kernel digests raw clouds directly
  // (1.2 s on 1024 points in 10D).
  if (dimension > 8)
    return cloud;
  // Qhull may fail on degenerate high-dim grids: fall back to the raw
  // cloud, which the exact kernel digests directly (slower, still exact)
  try
  {
    return ConvexHullMesher().build(cloud).getVertices();
  }
  catch (const std::exception &)
  {
    LOGWARN(OSS() << "IntersectionMesher hull reduction failed, using raw cloud");
    return cloud;
  }
}

// Reset-proof key read: ResourceMap::Reset() (run e.g. by Sphinx plot
// pre-code before every figure) wipes module runtime keys, whose static
// initializers run only once at load. Fall back to the registered default
// instead of throwing; core keys are unaffected (restored from defaults).
static UnsignedInteger GetUIntKey(const char * key,
                                  const UnsignedInteger fallback)
{
  try
  {
    return ResourceMap::GetAsUnsignedInteger(key);
  }
  catch (const Exception &)
  {
    return fallback;
  }
}

// Forward declaration (defined after buildConvex)
static Mesh AssembleConvexPieces(const Collection<Sample> & pieces,
                                 const UnsignedInteger dimension);

struct IntersectionMesherConvexTreePolicy
{
  const IntersectionMesher & intersectionMesher_;
  const Collection<Sample> & input_;
  Collection<Sample> & output_;

  IntersectionMesherConvexTreePolicy(const IntersectionMesher & intersectionMesher,
                                     const Collection<Sample> & input,
                                     Collection<Sample> & output)
    : intersectionMesher_(intersectionMesher)
    , input_(input)
    , output_(output)
  {}

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger n = r.begin(); n != r.end(); ++n)
      output_[n] = intersectionMesher_.buildConvexSample({input_[2 * n], input_[2 * n + 1]});
  }
};

// Support point of a vertex cloud along a direction (exact arithmetic)
static Point SupportPoint(const Sample & cloud,
                          const Point & direction)
{
  const UnsignedInteger size = cloud.getSize();
  UnsignedInteger best = 0;
  Scalar bestValue = cloud(0, 0) * direction[0];
  for (UnsignedInteger k = 1; k < cloud.getDimension(); ++ k)
    bestValue += cloud(0, k) * direction[k];
  for (UnsignedInteger i = 1; i < size; ++ i)
  {
    Scalar value = cloud(i, 0) * direction[0];
    for (UnsignedInteger k = 1; k < cloud.getDimension(); ++ k)
      value += cloud(i, k) * direction[k];
    if (value > bestValue)
    {
      bestValue = value;
      best = i;
    }
  }
  return cloud[best];
}

// Conservative disjointness witness for two convex clouds: a direction d
// with max(A.d) < min(B.d) proves disjointness (strict separation, exact
// floating-point support evaluations with tolerance margin). Never wrongly
// prunes: unproven pairs fall through to the exact kernel. Heuristic
// descent from the centroid difference, bounded iterations.
static Bool FindSeparatingDirection(const Sample & cloudA,
                                    const Sample & cloudB)
{
  const UnsignedInteger dimension = cloudA.getDimension();
  const Scalar scale = (cloudA.computeRange().norm() + cloudB.computeRange().norm()) / 2.0;
  const Scalar tolerance = std::sqrt(SpecFunc::Precision) * (1.0 + scale);
  Point direction(dimension);
  {
    const Point centerA(cloudA.computeMean());
    const Point centerB(cloudB.computeMean());
    Scalar norm = 0.0;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      direction[k] = centerA[k] - centerB[k];
      norm += direction[k] * direction[k];
    }
    if (!(norm > 0.0))
      return false;
    norm = std::sqrt(norm);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      direction[k] /= norm;
  }
  for (UnsignedInteger iteration = 0; iteration < 32; ++ iteration)
  {
    const Point supportA(SupportPoint(cloudA, direction));
    Point opposite(dimension);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      opposite[k] = -direction[k];
    const Point supportB(SupportPoint(cloudB, opposite));
    // gap along d: supportA.d - supportB.d < 0 proves separation
    Scalar gap = 0.0;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      gap += (supportA[k] - supportB[k]) * direction[k];
    if (gap < -tolerance)
      return true;
    // subgradient descent on the sphere
    Scalar gradNorm = 0.0;
    Point step(dimension);
    Scalar alignment = 0.0;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      step[k] = supportA[k] - supportB[k];
      alignment += step[k] * direction[k];
    }
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      step[k] -= alignment * direction[k];
      gradNorm += step[k] * step[k];
    }
    if (!(gradNorm > 0.0))
      return false;
    gradNorm = std::sqrt(gradNorm);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      direction[k] -= (0.5 / gradNorm) * step[k];
    Scalar norm = 0.0;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      norm += direction[k] * direction[k];
    if (!(norm > 0.0))
      return false;
    norm = std::sqrt(norm);
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      direction[k] /= norm;
  }
  return false;
}

// Bounding boxes of non-empty vertex clouds (caller-filtered)
static void ComputeCloudBounds(const Collection<Sample> & pieces,
                               Sample & lower,
                               Sample & upper)
{
  const UnsignedInteger size = pieces.getSize();
  const UnsignedInteger dimension = pieces[0].getDimension();
  lower = Sample(size, dimension);
  upper = Sample(size, dimension);
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    lower[i] = pieces[i].getMin();
    upper[i] = pieces[i].getMax();
  }
}

// Candidate (i, j) pairs with overlapping bounding boxes, as Indices of
// size 2. Exact filter: only truly disjoint pairs are skipped, so the
// downstream exact kernel sees the same non-empty pairs as the dense loop.
// Pieces spanning many cells are paired globally to bound grid insertions.
static Collection<Indices> ComputeCandidatePairs(const Collection<Sample> & pieces1,
                                                 const Collection<Sample> & pieces2)
{
  Collection<Indices> pairs(0);
  const UnsignedInteger n1 = pieces1.getSize();
  const UnsignedInteger n2 = pieces2.getSize();
  if (!n1 || !n2)
    return pairs;
  const UnsignedInteger dimension = pieces1[0].getDimension();
  if (!dimension)
  {
    for (UnsignedInteger i = 0; i < n1; ++ i)
      for (UnsignedInteger j = 0; j < n2; ++ j)
      {
        Indices pair(2);
        pair[0] = i;
        pair[1] = j;
        pairs.add(pair);
      }
    return pairs;
  }
  Sample lower1(0, dimension);
  Sample upper1(0, dimension);
  Sample lower2(0, dimension);
  Sample upper2(0, dimension);
  ComputeCloudBounds(pieces1, lower1, upper1);
  ComputeCloudBounds(pieces2, lower2, upper2);
  Point globalLower(dimension, SpecFunc::Infinity);
  Point globalUpper(dimension, -SpecFunc::Infinity);
  for (UnsignedInteger i = 0; i < n1; ++ i)
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      globalLower[k] = std::min(globalLower[k], lower1(i, k));
      globalUpper[k] = std::max(globalUpper[k], upper1(i, k));
    }
  for (UnsignedInteger j = 0; j < n2; ++ j)
    for (UnsignedInteger k = 0; k < dimension; ++ k)
    {
      globalLower[k] = std::min(globalLower[k], lower2(j, k));
      globalUpper[k] = std::max(globalUpper[k], upper2(j, k));
    }
  // cells per dimension for ~sqrt(n1*n2) cells total, capped for high dim
  UnsignedInteger cellsPerDim = static_cast<UnsignedInteger>(std::ceil(std::pow(1.0 * n1 * n2, 1.0 / (2 * dimension))));
  cellsPerDim = std::max(cellsPerDim, static_cast<UnsignedInteger>(1));
  UnsignedInteger totalCells = 1;
  for (UnsignedInteger k = 0; k < dimension; ++ k)
    totalCells *= cellsPerDim;
  while ((totalCells > 65536) && (cellsPerDim > 1))
  {
    cellsPerDim /= 2;
    totalCells = 1;
    for (UnsignedInteger k = 0; k < dimension; ++ k)
      totalCells *= cellsPerDim;
  }
  Point cellSize(dimension);
  for (UnsignedInteger k = 0; k < dimension; ++ k)
  {
    cellSize[k] = (globalUpper[k] - globalLower[k]) / cellsPerDim;
    if (!(cellSize[k] > 0.0))
      cellSize[k] = 1.0;
  }
  std::vector<UnsignedInteger> strides(dimension, 1);
  for (UnsignedInteger k = 1; k < dimension; ++ k)
    strides[k] = strides[k - 1] * cellsPerDim;
  std::vector< std::vector<uint64_t> > grid(totalCells);
  std::vector<UnsignedInteger> global1(0);
  std::vector<UnsignedInteger> global2(0);
  // encode (side, index) with side in the top bit region (indices < 2^56)
  for (UnsignedInteger pass = 0; pass < 2; ++ pass)
  {
    const Sample & lower = pass ? lower2 : lower1;
    const Sample & upper = pass ? upper2 : upper1;
    const UnsignedInteger count = pass ? n2 : n1;
    std::vector<UnsignedInteger> & globals = pass ? global2 : global1;
    for (UnsignedInteger i = 0; i < count; ++ i)
    {
      UnsignedInteger span = 1;
      std::vector<UnsignedInteger> first(dimension, 0);
      std::vector<UnsignedInteger> last(dimension, 0);
      for (UnsignedInteger k = 0; k < dimension; ++ k)
      {
        UnsignedInteger c0 = static_cast<UnsignedInteger>((lower(i, k) - globalLower[k]) / cellSize[k]);
        UnsignedInteger c1 = static_cast<UnsignedInteger>((upper(i, k) - globalLower[k]) / cellSize[k]);
        c0 = std::min(c0, cellsPerDim - 1);
        c1 = std::min(c1, cellsPerDim - 1);
        first[k] = c0;
        last[k] = c1;
        span *= (c1 - c0 + 1);
      }
      if (span > 64)
      {
        globals.push_back(i);
        continue;
      }
      // odometer over the spanned cell range
      std::vector<UnsignedInteger> cursor = first;
      for (;;)
      {
        UnsignedInteger cell = 0;
        for (UnsignedInteger k = 0; k < dimension; ++ k)
          cell += cursor[k] * strides[k];
        grid[cell].push_back((static_cast<uint64_t>(pass) << 56) | i);
        UnsignedInteger k = 0;
        while ((k < dimension) && (++cursor[k] > last[k]))
        {
          cursor[k] = first[k];
          ++ k;
        }
        if (k >= dimension)
          break;
      }
    }
  }
  std::vector<uint64_t> keys(0);
  for (UnsignedInteger g = 0; g < global1.size(); ++ g)
    for (UnsignedInteger j = 0; j < n2; ++ j)
      keys.push_back((static_cast<uint64_t>(global1[g]) << 32) | j);
  for (UnsignedInteger g = 0; g < global2.size(); ++ g)
    for (UnsignedInteger i = 0; i < n1; ++ i)
      keys.push_back((static_cast<uint64_t>(i) << 32) | global2[g]);
  for (UnsignedInteger cell = 0; cell < totalCells; ++ cell)
  {
    const std::vector<uint64_t> & content = grid[cell];
    for (UnsignedInteger a = 0; a < content.size(); ++ a)
    {
      if (content[a] >> 56)
        continue;
      const UnsignedInteger i = static_cast<UnsignedInteger>(content[a] & 0x00ffffffffffffffULL);
      for (UnsignedInteger b = 0; b < content.size(); ++ b)
      {
        if (!(content[b] >> 56))
          continue;
        const UnsignedInteger j = static_cast<UnsignedInteger>(content[b] & 0x00ffffffffffffffULL);
        // exact overlap re-test: sharing a cell is not sufficient
        Bool overlap = true;
        for (UnsignedInteger k = 0; k < dimension; ++ k)
          if ((upper1(i, k) < lower2(j, k)) || (upper2(j, k) < lower1(i, k)))
          {
            overlap = false;
            break;
          }
        if (overlap)
          keys.push_back((static_cast<uint64_t>(i) << 32) | j);
      }
    }
  }
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  for (UnsignedInteger q = 0; q < keys.size(); ++ q)
  {
    Indices pair(2);
    pair[0] = static_cast<UnsignedInteger>(keys[q] >> 32);
    pair[1] = static_cast<UnsignedInteger>(keys[q] & 0xffffffffu);
    pairs.add(pair);
  }
  return pairs;
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

struct IntersectionMesherIndexedPolicy
{
  const IntersectionMesher & intersectionMesher_;
  const Collection<Sample> & input1_;
  const Collection<Sample> & input2_;
  const Collection<Indices> & pairs_;
  Collection<Sample> & output_;

  IntersectionMesherIndexedPolicy(const IntersectionMesher & intersectionMesher,
                                  const Collection<Sample> & input1,
                                  const Collection<Sample> & input2,
                                  const Collection<Indices> & pairs,
                                  Collection<Sample> & output)
    : intersectionMesher_(intersectionMesher)
    , input1_(input1)
    , input2_(input2)
    , pairs_(pairs)
    , output_(output)
  {}

  inline void operator()(const TBBImplementation::BlockedRange<UnsignedInteger> & r) const
  {
    for (UnsignedInteger n = r.begin(); n != r.end(); ++n)
      output_[n] = intersectionMesher_.buildConvexSample({input1_[pairs_[n][0]], input2_[pairs_[n][1]]});
  }
};

// Pairwise intersection of two convex-piece lists: dense blocked loop for
// small products, AABB-grid candidate pairs above the GridThreshold.
// Returns the non-empty intersections.
static Collection<Sample> IntersectPieceLists(const IntersectionMesher & mesher,
                                              const Collection<Sample> & pieces1,
                                              const Collection<Sample> & pieces2)
{
  Collection<Sample> result(0);
  // drop empty clouds (union members contributing nothing); an emptied list
  // empties the pairwise product
  Collection<Sample> dense1(0);
  for (UnsignedInteger i = 0; i < pieces1.getSize(); ++ i)
    if (pieces1[i].getSize())
      dense1.add(pieces1[i]);
  Collection<Sample> dense2(0);
  for (UnsignedInteger j = 0; j < pieces2.getSize(); ++ j)
    if (pieces2[j].getSize())
      dense2.add(pieces2[j]);
  if (!dense1.getSize() || !dense2.getSize())
    return result;
  const Collection<Sample> & kept1 = dense1;
  const Collection<Sample> & kept2 = dense2;
  const UnsignedInteger toDoSize = kept1.getSize() * kept2.getSize();
  const UnsignedInteger blockSize = GetUIntKey("IntersectionMesher-BlockSize", 1 << 16);
  const UnsignedInteger gridThreshold = GetUIntKey("IntersectionMesher-GridThreshold", 1 << 10);
  if (toDoSize > gridThreshold)
  {
    const Collection<Indices> pairs(ComputeCandidatePairs(kept1, kept2));
    const UnsignedInteger pairsSize = pairs.getSize();
    LOGDEBUG(OSS() << "Grid filter: " << toDoSize << " pairs -> " << pairsSize << " candidates");
    Collection<Sample> resultChunk(pairsSize);
    for (UnsignedInteger done = 0; done < pairsSize; done += blockSize)
    {
      const UnsignedInteger actualBlockSize = std::min(blockSize, pairsSize - done);
      const IntersectionMesherIndexedPolicy policy(mesher, kept1, kept2, pairs, resultChunk);
      TBBImplementation::ParallelFor(done, done + actualBlockSize, policy);
      // prune empty intersections
      for (UnsignedInteger i0 = done; i0 < done + actualBlockSize; ++ i0)
        if (resultChunk[i0].getSize())
          result.add(resultChunk[i0]);
    }
    return result;
  }
  Collection<Sample> denseChunk(blockSize);
  for (UnsignedInteger done = 0; done < toDoSize; done += blockSize)
  {
    const UnsignedInteger actualBlockSize = std::min(blockSize, toDoSize - done);
    const IntersectionMesherConvexSamplePolicy policy(mesher, kept1, kept2, done, denseChunk);
    TBBImplementation::ParallelFor(0, actualBlockSize, policy);
    // prune empty intersections
    for (UnsignedInteger i0 = 0; i0 < actualBlockSize; ++ i0)
      if (denseChunk[i0].getSize())
        result.add(denseChunk[i0]);
  }
  return result;
}

Mesh IntersectionMesher::build(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

  // build decomposition of first mesh. Convex inputs are kept whole (one
  // convex piece) instead of split into simplices: fewer pairs downstream.
  ConvexDecompositionMesher convexDecompositionMesher;
  convexDecompositionMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
#if 0
  LOGTRACE(OSS() << "Build decomposition of mesh 0");
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
#endif
  Collection<Sample> unionCurrent;
  if (coll[0].getSimplicesNumber() && (coll[0].isConvex() || ConvexDecompositionMesher::IsConvex(coll[0])))
    unionCurrent.add(ReduceConvexCloud(coll[0].getVertices()));
  else
  {
    const Collection<Mesh> baseDecomposition0(convexDecompositionMesher.build(coll[0]));
    const UnsignedInteger baseDecompositionSize0 = baseDecomposition0.getSize();
    for (UnsignedInteger j = 0; j < baseDecompositionSize0; ++ j)
      unionCurrent.add(baseDecomposition0[j].getVertices());
  }
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
    if (coll[i].getSimplicesNumber() && (coll[i].isConvex() || ConvexDecompositionMesher::IsConvex(coll[i])))
      unionNext.add(ReduceConvexCloud(coll[i].getVertices()));
    else
    {
      const Collection<Mesh> baseDecomposition(convexDecompositionMesher.build(coll[i]));
      const UnsignedInteger baseDecompositionSize = baseDecomposition.getSize();
      //LOGTRACE(OSS() << "Got " << baseDecompositionSize << " parts");
      for (UnsignedInteger j = 0; j < baseDecompositionSize; ++ j)
        unionNext.add(baseDecomposition[j].getVertices());
    }
#if 0
    t1 = std::chrono::steady_clock::now();
    timeDuration = std::chrono::duration<Scalar>(t1 - t0).count();
    LOGTRACE(OSS() << "Got " << baseDecompositionSize << " parts in t=" << timeDuration << "s");
    t0 = std::chrono::steady_clock::now();
#endif
    // loop over intersections (dense or AABB-grid candidates)
    const Collection<Sample> result(IntersectPieceLists(*this, unionCurrent, unionNext));
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

  // assemble the surviving convex pieces (single piece: hull+volume mesher)
  return AssembleConvexPieces(unionCurrent, coll[0].getDimension());
}

Mesh IntersectionMesher::buildWithConvexParts(const Mesh & mesh, const SampleCollection & convexPieces) const
{
  ConvexDecompositionMesher convexDecompositionMesher;
  convexDecompositionMesher.setUseSimplicesDecomposition(useSimplicesDecomposition_);
  Collection<Sample> unionCurrent;
  const Collection<Mesh> baseDecomposition0(convexDecompositionMesher.build(mesh));
  for (UnsignedInteger j = 0; j < baseDecomposition0.getSize(); ++ j)
    unionCurrent.add(baseDecomposition0[j].getVertices());

  const Collection<Sample> result(IntersectPieceLists(*this, unionCurrent, convexPieces));
  return AssembleConvexPieces(result, mesh.getDimension());
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


// Assemble convex pieces into a mesh: single piece via hull+volume mesher
// (few tetrahedra), several via per-piece Delaunay plus union
static Mesh AssembleConvexPieces(const Collection<Sample> & pieces,
                                 const UnsignedInteger dimension)
{
  if (!pieces.getSize())
    return Mesh(Sample(0, dimension));
  if (pieces.getSize() == 1)
    return VolumeMesher().build(ConvexHullMesher().build(pieces[0]));
  CloudMesher cloudMesher;
  Collection<Mesh> collMesh(pieces.getSize());
  for (UnsignedInteger i = 0; i < pieces.getSize(); ++ i)
    collMesh[i] = cloudMesher.build(pieces[i]);
  return UnionMesher().build(collMesh);
}

Mesh IntersectionMesher::buildConvex(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];
  // pieces per mesh (union semantics within a mesh), folded pairwise across
  // meshes: split pieces are disjuncts and must flow through the cartesian
  // product plus union, never through a conjunctive tree/mono reduction.
  // An empty piece list (empty mesh) empties the whole intersection.
  const UnsignedInteger dimension = coll[0].getDimension();
  const Sample firstPiece(ReduceConvexCloud(coll[0].getVertices()));
  Collection<Sample> current(0);
  if (firstPiece.getSize())
    current.add(firstPiece);
  if (!current.getSize())
    return Mesh(Sample(0, dimension));
  for (UnsignedInteger i = 1; i < size; ++ i)
  {
    const Sample nextPiece(ReduceConvexCloud(coll[i].getVertices()));
    Collection<Sample> next(0);
    if (nextPiece.getSize())
      next.add(nextPiece);
    if (!next.getSize())
      return Mesh(Sample(0, dimension));
    current = IntersectPieceLists(*this, current, next);
    if (!current.getSize())
      return Mesh(Sample(0, dimension));
  }
  return AssembleConvexPieces(current, dimension);
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
    // A pruned input is disjoint from (or only touching) the running box,
    // hence from the running intersection: the total intersection is empty
    // (as a volume; boundary-only contacts carry no d-volume to mesh).
    // Returning now also avoids enumerating the survivors, which would
    // wrongly produce their non-empty intersection.
    if (pruned)
      return Sample(0, dimension);
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
  // Conservative exact-disjointness filter: pairwise separating directions
  // skip the exact enumeration for well-separated clouds (never prunes
  // wrongly; unproven pairs fall through). Runs on bbox survivors only.
  for (UnsignedInteger a = 0; a < remainingSize; ++ a)
    for (UnsignedInteger b = a + 1; b < remainingSize; ++ b)
      if (FindSeparatingDirection(coll[remainingIndices[a]], coll[remainingIndices[b]]))
        return Sample(0, dimension);

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


// Binary-tree pairwise reduction of convex pieces: pairs intersect in
// parallel per level, odd piece carried over, empties pruned per level
// with early exit. Each enumeration stays small while intermediates
// shrink, unlike the monolithic single enumeration of buildConvexSample.
Sample IntersectionMesher::buildConvexTree(const SampleCollection & coll) const
{
  const UnsignedInteger dimension = coll.getSize() ? coll[0].getDimension() : 0;
  // NB: inputs are conjuncts (all intersected); only same-polytope hull
  // reduction applies here, never splitting (split pieces are disjuncts)
  Collection<Sample> current(0);
  for (UnsignedInteger i = 0; i < coll.getSize(); ++ i)
    if (coll[i].getSize())
      current.add(ReduceConvexCloud(coll[i]));
  if (current.getSize() < 2)
    return current.getSize() ? current[0] : Sample(0, dimension);
  while (current.getSize() > 1)
  {
    const UnsignedInteger pairsNumber = current.getSize() / 2;
    const UnsignedInteger oddNumber = current.getSize() % 2;
    Collection<Sample> next(pairsNumber + oddNumber);
    const IntersectionMesherConvexTreePolicy policy(*this, current, next);
    TBBImplementation::ParallelFor(0, pairsNumber, policy);
    if (oddNumber)
      next[pairsNumber] = current[current.getSize() - 1];
    current.resize(0);
    for (UnsignedInteger i = 0; i < next.getSize(); ++ i)
      if (next[i].getSize())
        current.add(next[i]);
    if (!current.getSize())
      return Sample(0, dimension);
  }
  return current[0];
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

    // loop over intersections (dense or AABB-grid candidates)
    unionCurrent = IntersectPieceLists(*this, unionCurrent, unionNext);

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
    ResourceMap::AddAsUnsignedInteger("IntersectionMesher-GridThreshold", 1 << 10);
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
