//                                               -*- C++ -*-
/**
 *  @brief Union meshing
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
#include "otmeshing/UnionMesher.hxx"

#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/SpecFunc.hxx>
#include <openturns/KDTree.hxx>

#include <algorithm>
#include <vector>

using namespace OT;

namespace OTMESHING
{

CLASSNAMEINIT(UnionMesher)
static const Factory<UnionMesher> Factory_UnionMesher;


/* Default constructor */
UnionMesher::UnionMesher()
  : PersistentObject()
{
  // Nothing to do
}

/* Virtual constructor */
UnionMesher * UnionMesher::clone() const
{
  return new UnionMesher(*this);
}

/* String converter */
String UnionMesher::__repr__() const
{
  OSS oss(true);
  oss << "class=" << UnionMesher::GetClassName();
  return oss;
}

/* String converter */
String UnionMesher::__str__(const String & ) const
{
  return __repr__();
}

/* Deduplicate mesh vertices */
Mesh UnionMesher::CompressMesh(const Mesh & mesh)
{
  const UnsignedInteger dimension = mesh.getDimension();
  const Sample vertices(mesh.getVertices());
  const UnsignedInteger fullSize = vertices.getSize();
  if (!fullSize)
    return mesh;
  IndicesCollection simplices(mesh.getSimplices());

  // mark vertices referenced by a simplex
  Indices usedVertex(fullSize, 0);
  for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    for (UnsignedInteger j = 0; j <= dimension; ++ j)
      usedVertex[simplices(i, j)] = 1;

  // Phase 1: union-find to build connected components of vertices within tolerance.
  // This ensures transitive closure: if A is close to B and B is close to C,
  // all three end up in the same component even if A and C are not directly close.
  Indices parent(fullSize);
  parent.fill();

  // iterative find with full path compression
  auto find = [&](UnsignedInteger x) -> UnsignedInteger
  {
    UnsignedInteger root = x;
    while (root != parent[root])
      root = parent[root];
    while (x != root)
    {
      const UnsignedInteger next = parent[x];
      parent[x] = root;
      x = next;
    }
    return root;
  };

  const KDTree tree(vertices);
  const Scalar tolerance = SpecFunc::Precision * vertices.computeRange().norm();
  for (UnsignedInteger i = 0; i < fullSize; ++ i)
  {
    if (!usedVertex[i])
      continue;
    Point distance;
    const Indices nearest(tree.queryRadius(vertices[i], tolerance, distance));
    for (UnsignedInteger k = 0; k < nearest.getSize(); ++ k)
    {
      const UnsignedInteger j = nearest[k];
      if (!usedVertex[j] || j == i)
        continue;
      // Recompute rootI on each iteration so transitive chains merge correctly
      const UnsignedInteger rootI = find(i);
      const UnsignedInteger rootJ = find(j);
      if (rootI != rootJ)
        parent[rootI] = rootJ;
    }
  }

  // Phase 2: identify unique roots and assign compressed indices
  Indices compressedVertexMap(fullSize, fullSize);
  UnsignedInteger nRoots = 0;
  for (UnsignedInteger i = 0; i < fullSize; ++ i)
  {
    if (!usedVertex[i])
      continue;
    const UnsignedInteger r = find(i);
    if (compressedVertexMap[r] >= fullSize)
    {
      compressedVertexMap[r] = nRoots;
      ++ nRoots;
    }
  }

  // Phase 3: accumulate component sums into the final sample
  Sample verticesCompressed(nRoots, dimension);
  Indices sizes(nRoots, 0);
  for (UnsignedInteger i = 0; i < fullSize; ++ i)
  {
    if (!usedVertex[i])
      continue;
    const UnsignedInteger r = find(i);
    const UnsignedInteger idx = compressedVertexMap[r];
    for (UnsignedInteger d = 0; d < dimension; ++ d)
      verticesCompressed[idx][d] += vertices[i][d];
    sizes[idx] += 1;
  }

  // Phase 4: compute centroid of each component
  for (UnsignedInteger i = 0; i < nRoots; ++ i)
  {
    const Scalar invSize = 1.0 / sizes[i];
    for (UnsignedInteger d = 0; d < dimension; ++ d)
      verticesCompressed[i][d] *= invSize;
  }

  // Phase 5: build full mapping from original index to compressed index
  for (UnsignedInteger i = 0; i < fullSize; ++ i)
  {
    if (!usedVertex[i])
      continue;
    const UnsignedInteger r = find(i);
    compressedVertexMap[i] = compressedVertexMap[r];
  }

  LOGDEBUG(OSS() << "recompression fullSize=" << fullSize << " compressedSize=" << nRoots);

  // renumber vertex indices
  for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    for (UnsignedInteger j = 0; j <= dimension; ++ j)
      simplices(i, j) = compressedVertexMap[simplices(i, j)];
  return Mesh(verticesCompressed, simplices);
}

namespace
{
// Reset-proof key read: see IntersectionMesher.cxx (Sphinx plot pre-code
// calls ResourceMap::Reset(), wiping module runtime keys).
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

// Sorted vertex key of a simplex, for exact-duplicate detection
Indices SortedKey(const IndicesCollection & simplices,
                  const UnsignedInteger simplexIndex,
                  const UnsignedInteger verticesPerSimplex)
{
  Indices key(verticesPerSimplex);
  for (UnsignedInteger j = 0; j < verticesPerSimplex; ++ j)
    key[j] = simplices(simplexIndex, j);
  std::sort(key.begin(), key.end());
  return key;
}

// Strictly-inside point-on-segment test in any dimension.
static Bool SplitPointOnEdgeND(const Point & p,
                               const Point & a,
                               const Point & b,
                               const Scalar tolerance,
                               Scalar & splitParameter)
{
  const UnsignedInteger dimension = p.getDimension();
  Scalar squaredLength = 0.0;
  Scalar dot = 0.0;
  for (UnsignedInteger k = 0; k < dimension; ++ k)
  {
    const Scalar ab = b[k] - a[k];
    squaredLength += ab * ab;
    dot += (p[k] - a[k]) * ab;
  }
  if (!(squaredLength > 0.0))
    return false;
  const Scalar t = dot / squaredLength;
  const Scalar edgeLength = std::sqrt(squaredLength);
  if (!(t * edgeLength > tolerance) || !((1.0 - t) * edgeLength > tolerance))
    return false;
  Scalar squaredDistance = 0.0;
  for (UnsignedInteger k = 0; k < dimension; ++ k)
  {
    const Scalar delta = p[k] - (a[k] + t * (b[k] - a[k]));
    squaredDistance += delta * delta;
  }
  if (!(std::sqrt(squaredDistance) <= tolerance))
    return false;
  splitParameter = t;
  return true;
}

// Strictly-inside point-on-triangle test in 3D (face T-junction detection).
// Barycentric coordinates on the dominant plane, strict margins.
static Bool SplitPointOnTriangle(const Point & p,
                                 const Point & a,
                                 const Point & b,
                                 const Point & c,
                                 const Scalar tolerance)
{
  // dominant axis of the facet normal
  const Scalar ux = b[0] - a[0];
  const Scalar uy = b[1] - a[1];
  const Scalar uz = b[2] - a[2];
  const Scalar vx = c[0] - a[0];
  const Scalar vy = c[1] - a[1];
  const Scalar vz = c[2] - a[2];
  const Scalar nx = uy * vz - uz * vy;
  const Scalar ny = uz * vx - ux * vz;
  const Scalar nz = ux * vy - uy * vx;
  const Scalar anx = std::abs(nx);
  const Scalar any = std::abs(ny);
  const Scalar anz = std::abs(nz);
  UnsignedInteger i = 0;
  UnsignedInteger j = 1;
  if ((anx >= any) && (anx >= anz))
  {
    i = 1;
    j = 2;
  }
  else if ((any >= anx) && (any >= anz))
  {
    i = 0;
    j = 2;
  }
  // 2D solve on axes (i, j)
  const Scalar d00 = (b[i] - a[i]) * (b[i] - a[i]) + (b[j] - a[j]) * (b[j] - a[j]);
  const Scalar d01 = (b[i] - a[i]) * (c[i] - a[i]) + (b[j] - a[j]) * (c[j] - a[j]);
  const Scalar d11 = (c[i] - a[i]) * (c[i] - a[i]) + (c[j] - a[j]) * (c[j] - a[j]);
  const Scalar det = d00 * d11 - d01 * d01;
  if (!(std::abs(det) > 0.0))
    return false;
  const Scalar d20 = (p[i] - a[i]) * (b[i] - a[i]) + (p[j] - a[j]) * (b[j] - a[j]);
  const Scalar d21 = (p[i] - a[i]) * (c[i] - a[i]) + (p[j] - a[j]) * (c[j] - a[j]);
  const Scalar v = (d11 * d20 - d01 * d21) / det;
  const Scalar w = (d00 * d21 - d01 * d20) / det;
  const Scalar u = 1.0 - v - w;
  const Scalar margin = std::sqrt(SpecFunc::Precision);
  if (!(u > margin) || !(v > margin) || !(w > margin))
    return false;
  // distance to the facet plane
  const Scalar norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (!(norm > 0.0))
    return false;
  const Scalar distance = std::abs(nx * (p[0] - a[0]) + ny * (p[1] - a[1]) + nz * (p[2] - a[2])) / norm;
  return distance <= tolerance;
}

// Strictly-inside point-on-segment test in 2D (T-junction detection).
// Returns true and sets splitParameter when point p lies on the open
// segment (a, b) within tolerance.
Bool SplitPointOnEdge(const Point & p,
                      const Point & a,
                      const Point & b,
                      const Scalar tolerance,
                      Scalar & splitParameter)
{
  const Scalar abx = b[0] - a[0];
  const Scalar aby = b[1] - a[1];
  const Scalar squaredLength = abx * abx + aby * aby;
  if (!(squaredLength > 0.0))
    return false;
  const Scalar t = ((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / squaredLength;
  const Scalar edgeLength = std::sqrt(squaredLength);
  if (!(t * edgeLength > tolerance) || !((1.0 - t) * edgeLength > tolerance))
    return false;
  const Scalar dx = p[0] - (a[0] + t * abx);
  const Scalar dy = p[1] - (a[1] + t * aby);
  if (!(std::sqrt(dx * dx + dy * dy) <= tolerance))
    return false;
  splitParameter = t;
  return true;
}
} // anonymous namespace

Mesh UnionMesher::build(const MeshCollection & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return CompressMesh(coll[0]);

  const UnsignedInteger dimension = coll[0].getDimension();
  for (UnsignedInteger i = 1; i < size; ++ i)
    if (coll[i].getDimension() != dimension)
      throw InvalidArgumentException(HERE) << "UnionMesher expected meshes of same dimension";

  UnsignedInteger simplicesNumber = 0;
  Sample vertices(0, dimension);
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    vertices.add(coll[i].getVertices());
    simplicesNumber += coll[i].getSimplicesNumber();
  }
  IndicesCollection simplices(simplicesNumber, dimension + 1);
  UnsignedInteger simplicesOffset = 0;
  UnsignedInteger vertexOffset = 0;
  for (UnsignedInteger i = 0; i < size; ++ i)
  {
    IndicesCollection simplicesI(coll[i].getSimplices());
    const UnsignedInteger sizeI = simplicesI.getSize();
    for (UnsignedInteger j = 0; j < sizeI; ++ j)
      for (UnsignedInteger k = 0; k <= dimension; ++ k)
        simplices(simplicesOffset + j, k) = simplicesI(j, k) + vertexOffset;
    simplicesOffset += sizeI;
    vertexOffset += coll[i].getVerticesNumber();
  }
  const Mesh welded(CompressMesh(Mesh(vertices, simplices, false)));
  const UnsignedInteger weldedSimplicesNumber = welded.getSimplicesNumber();
  if (!weldedSimplicesNumber)
    return welded;
  const Sample weldedVertices(welded.getVertices());
  const IndicesCollection weldedSimplices(welded.getSimplices());
  const UnsignedInteger verticesPerSimplex = dimension + 1;
  // U1: drop exact-duplicate simplices (same welded vertex set), keeping the
  // first occurrence. Coincident inputs then union to their common mesh
  // instead of double-counting volumes. Sort-based for O(n log n).
  std::vector< std::pair<Indices, UnsignedInteger> > keyedSimplices(weldedSimplicesNumber);
  for (UnsignedInteger i = 0; i < weldedSimplicesNumber; ++ i)
  {
    keyedSimplices[i].first = SortedKey(weldedSimplices, i, verticesPerSimplex);
    keyedSimplices[i].second = i;
  }
  std::sort(keyedSimplices.begin(), keyedSimplices.end());
  Indices keptSimplices(0);
  for (UnsignedInteger u = 0; u < keyedSimplices.size(); ++ u)
    if ((u == 0) || (keyedSimplices[u].first != keyedSimplices[u - 1].first))
      keptSimplices.add(keyedSimplices[u].second);
  LOGDEBUG(OSS() << "UnionMesher dropped " << (weldedSimplicesNumber - keptSimplices.getSize()) << " duplicate simplices");
  // U2 (2D): resolve T-junctions by subdivision. A welded vertex strictly
  // inside an edge of a triangle it does not belong to (meshes sharing
  // vertices connected differently) splits that triangle, so interfaces
  // match node-for-node instead of keeping hanging nodes.
  Collection<Indices> finalSimplices(0);
  for (UnsignedInteger s = 0; s < keptSimplices.getSize(); ++ s)
  {
    Indices simplex(verticesPerSimplex);
    for (UnsignedInteger j = 0; j < verticesPerSimplex; ++ j)
      simplex[j] = weldedSimplices(keptSimplices[s], j);
    finalSimplices.add(simplex);
  }
  // T-junction refinement is an interface repair for few-piece unions. Many-piece
  // gluings (e.g. IntersectionMesher assembly of overlapping convex pieces) keep
  // weld plus duplicate removal only: subdivision would cascade on volumetric
  // overlaps instead of converging on interfaces.
  const UnsignedInteger subdivisionThreshold = GetUIntKey("UnionMesher-SubdivisionThreshold", 256);
  const Bool refineInterfaces = (finalSimplices.getSize() <= subdivisionThreshold);
  LOGDEBUG(OSS() << "UnionMesher pieces=" << finalSimplices.getSize() << " refine=" << refineInterfaces);
  if (refineInterfaces && (dimension == 2))
  {
    const Scalar tolerance = SpecFunc::Precision * weldedVertices.computeRange().norm();
    const UnsignedInteger verticesNumber = weldedVertices.getSize();
    // worklist of triangle vertex triplets, refined until no split applies
    UnsignedInteger guard = 0;
    const UnsignedInteger guardLimit = 4 * (verticesNumber + finalSimplices.getSize()) + 1;
    Bool refined = true;
    while (refined && (guard < guardLimit))
    {
      refined = false;
      ++ guard;
      for (UnsignedInteger t = 0; t < finalSimplices.getSize(); ++ t)
      {
        const Indices tri(finalSimplices[t]);
        // split points per edge: (edge position 0..2, vertex, parameter)
        UnsignedInteger splitEdge = 3;
        UnsignedInteger splitVertex = verticesNumber;
        for (UnsignedInteger e = 0; (e < 3) && (splitEdge > 2); ++ e)
        {
          const UnsignedInteger a = tri[e];
          const UnsignedInteger b = tri[(e + 1) % 3];
          const Point pointA(weldedVertices[a]);
          const Point pointB(weldedVertices[b]);
          for (UnsignedInteger m = 0; m < verticesNumber; ++ m)
          {
            if ((m == a) || (m == b))
              continue;
            Scalar unusedParameter = 0.0;
            if (SplitPointOnEdge(weldedVertices[m], pointA, pointB, tolerance, unusedParameter))
            {
              splitEdge = e;
              splitVertex = m;
              break;
            }
          }
        }
        if (splitEdge > 2)
          continue;
        // 1-to-2 split of tri along (a, b) at splitVertex
        const UnsignedInteger a = tri[splitEdge];
        const UnsignedInteger b = tri[(splitEdge + 1) % 3];
        const UnsignedInteger w = tri[(splitEdge + 2) % 3];
        Indices first(3);
        first[0] = a;
        first[1] = splitVertex;
        first[2] = w;
        Indices second(3);
        second[0] = splitVertex;
        second[1] = b;
        second[2] = w;
        finalSimplices[t] = first;
        finalSimplices.add(second);
        refined = true;
      }
    }
    LOGDEBUG(OSS() << "UnionMesher 2D subdivision passes=" << guard);
  }
  if (refineInterfaces && (dimension == 3))
  {
    const Scalar tolerance = SpecFunc::Precision * weldedVertices.computeRange().norm();
    const UnsignedInteger verticesNumber = weldedVertices.getSize();
    UnsignedInteger guard = 0;
    const UnsignedInteger guardLimit = 4 * (verticesNumber + finalSimplices.getSize()) + 1;
    Bool refined = true;
    while (refined && (guard < guardLimit))
    {
      refined = false;
      ++ guard;
      for (UnsignedInteger t = 0; t < finalSimplices.getSize(); ++ t)
      {
        const Indices tet(finalSimplices[t]);
        // edge splits first (1-to-2), then face splits (1-to-3)
        UnsignedInteger splitEdgeA = verticesNumber;
        UnsignedInteger splitEdgeB = verticesNumber;
        UnsignedInteger splitMid = verticesNumber;
        for (UnsignedInteger e1 = 0; (e1 < 4) && (splitMid >= verticesNumber); ++ e1)
          for (UnsignedInteger e2 = e1 + 1; (e2 < 4) && (splitMid >= verticesNumber); ++ e2)
          {
            const UnsignedInteger a = tet[e1];
            const UnsignedInteger b = tet[e2];
            for (UnsignedInteger m = 0; m < verticesNumber; ++ m)
            {
              if ((m == a) || (m == b))
                continue;
              Scalar unusedParameter = 0.0;
              if (SplitPointOnEdgeND(weldedVertices[m], weldedVertices[a], weldedVertices[b], tolerance, unusedParameter))
              {
                splitEdgeA = a;
                splitEdgeB = b;
                splitMid = m;
                break;
              }
            }
          }
        if (splitMid < verticesNumber)
        {
          // 1-to-2 split: replace tet by the two halves sharing (mid, apex0, apex1)
          UnsignedInteger apex[2] = {verticesNumber, verticesNumber};
          UnsignedInteger apexCount = 0;
          for (UnsignedInteger e = 0; e < 4; ++ e)
            if ((tet[e] != splitEdgeA) && (tet[e] != splitEdgeB))
              apex[apexCount++] = tet[e];
          Indices first(4);
          first[0] = splitEdgeA;
          first[1] = splitMid;
          first[2] = apex[0];
          first[3] = apex[1];
          Indices second(4);
          second[0] = splitMid;
          second[1] = splitEdgeB;
          second[2] = apex[0];
          second[3] = apex[1];
          finalSimplices[t] = first;
          finalSimplices.add(second);
          refined = true;
          continue;
        }
        // face splits (1-to-3 through the face point)
        UnsignedInteger splitFace[3] = {verticesNumber, verticesNumber, verticesNumber};
        UnsignedInteger faceMid = verticesNumber;
        for (UnsignedInteger f = 0; (f < 4) && (faceMid >= verticesNumber); ++ f)
        {
          // facet opposite tet[f]: the other three vertices
          UnsignedInteger fa = verticesNumber;
          UnsignedInteger fb = verticesNumber;
          UnsignedInteger fc = verticesNumber;
          UnsignedInteger seen = 0;
          for (UnsignedInteger e = 0; e < 4; ++ e)
            if (e != f)
            {
              if (seen == 0) fa = tet[e];
              else if (seen == 1) fb = tet[e];
              else fc = tet[e];
              ++ seen;
            }
          for (UnsignedInteger m = 0; m < verticesNumber; ++ m)
          {
            if ((m == fa) || (m == fb) || (m == fc))
              continue;
            if (SplitPointOnTriangle(weldedVertices[m], weldedVertices[fa], weldedVertices[fb], weldedVertices[fc], tolerance))
            {
              splitFace[0] = fa;
              splitFace[1] = fb;
              splitFace[2] = fc;
              faceMid = m;
              break;
            }
          }
        }
        if (faceMid < verticesNumber)
        {
          UnsignedInteger apex = verticesNumber;
          for (UnsignedInteger e = 0; e < 4; ++ e)
            if ((tet[e] != splitFace[0]) && (tet[e] != splitFace[1]) && (tet[e] != splitFace[2]))
              apex = tet[e];
          if (apex >= verticesNumber)
            continue;
          Indices first(4);
          first[0] = splitFace[0];
          first[1] = splitFace[1];
          first[2] = faceMid;
          first[3] = apex;
          Indices second(4);
          second[0] = splitFace[1];
          second[1] = splitFace[2];
          second[2] = faceMid;
          second[3] = apex;
          Indices third(4);
          third[0] = splitFace[2];
          third[1] = splitFace[0];
          third[2] = faceMid;
          third[3] = apex;
          finalSimplices[t] = first;
          finalSimplices.add(second);
          finalSimplices.add(third);
          refined = true;
        }
      }
    }
    LOGDEBUG(OSS() << "UnionMesher 3D subdivision passes=" << guard);
  }
  // compact unused vertices and return a checked mesh
  const UnsignedInteger finalSize = finalSimplices.getSize();
  Indices usedVertex(weldedVertices.getSize(), 0);
  for (UnsignedInteger i = 0; i < finalSize; ++ i)
    for (UnsignedInteger j = 0; j < verticesPerSimplex; ++ j)
      usedVertex[finalSimplices[i][j]] = 1;
  Indices compactMap(weldedVertices.getSize(), weldedVertices.getSize());
  UnsignedInteger compactSize = 0;
  for (UnsignedInteger i = 0; i < weldedVertices.getSize(); ++ i)
    if (usedVertex[i])
      compactMap[i] = compactSize++;
  Sample compactVertices(compactSize, dimension);
  for (UnsignedInteger i = 0; i < weldedVertices.getSize(); ++ i)
    if (usedVertex[i])
      compactVertices[compactMap[i]] = weldedVertices[i];
  IndicesCollection compactSimplices(finalSize, verticesPerSimplex);
  for (UnsignedInteger i = 0; i < finalSize; ++ i)
    for (UnsignedInteger j = 0; j < verticesPerSimplex; ++ j)
      compactSimplices(i, j) = compactMap[finalSimplices[i][j]];
  return Mesh(compactVertices, compactSimplices);
}

struct UnionMesher_init
{
  UnionMesher_init()
  {
    ResourceMap::AddAsUnsignedInteger("UnionMesher-SubdivisionThreshold", 256);
  }
};

static UnionMesher_init __UnionMesher_initializer;

/* Method save() stores the object through the StorageManager */
void UnionMesher::save(Advocate & adv) const
{
  PersistentObject::save(adv);
}

/* Method load() reloads the object from the StorageManager */
void UnionMesher::load(Advocate & adv)
{
  PersistentObject::load(adv);
}

}
