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
  return CompressMesh(Mesh(vertices, simplices, false));
}

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
