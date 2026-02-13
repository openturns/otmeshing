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
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/KDTree.hxx"

#include "otmeshing/IntersectionMesher.hxx"
#include "otmeshing/CloudMesher.hxx"

#ifdef OPENTURNS_HAVE_CDDLIB
#include <setoper.h>
#include <cdd.h>
#endif

#include <nanoflann.hpp>
#if NANOFLANN_VERSION < 0x150
namespace nanoflann
{
using SearchParameters = SearchParams;
}
#endif

using namespace OT;

namespace OTMESHING
{

class KDTreeSampleAdaptor
{
public:
  explicit KDTreeSampleAdaptor(const Sample & points)
    : data_(points.getImplementation()->data())
    , size_(points.getSize())
    , dimension_(points.getDimension())
  {
  }

  inline size_t kdtree_get_point_count() const
  {
    return size_;
  }

  inline Scalar kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    return data_[dim + idx * dimension_];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /*bb*/) const
  {
    return false;
  }

private:
  const Scalar *data_ = nullptr;
  UnsignedInteger size_ = 0;
  UnsignedInteger dimension_ = 0;
};

using nano_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor <
                       nanoflann::L2_Simple_Adaptor<Scalar, KDTreeSampleAdaptor >,
                       KDTreeSampleAdaptor, -1 >;

class KDTree2
{
public:
  explicit KDTree2(const Sample & points)
  {
    const UnsignedInteger dimension = points.getDimension();
    sampleAdaptor_ = new KDTreeSampleAdaptor(points);
    nanoflann::KDTreeSingleIndexAdaptorParams indexParameters;
    indexParameters.leaf_max_size = ResourceMap::GetAsUnsignedInteger("KDTree-leaf_max_size");
#if NANOFLANN_VERSION >= 0x150
    indexParameters.n_thread_build = ResourceMap::GetAsUnsignedInteger("KDTree-n_thread_build");
#endif
    indexAdaptor_ = new nano_kd_tree_t(dimension, *sampleAdaptor_, indexParameters);
  }

  Indices queryRadius(const Point & x, const Scalar radius, Point & distanceOut, const Bool sorted = false) const
  {
#if NANOFLANN_VERSION >= 0x150
    std::vector<nanoflann::ResultItem<unsigned int, Scalar> > indicesDists;
#else
    std::vector<std::pair<UnsignedInteger, Scalar> > indicesDists;
#endif
    nanoflann::SearchParameters searchParameters;
    searchParameters.sorted = sorted;
    const UnsignedInteger nFound = indexAdaptor_->radiusSearch(x.data(), radius * radius, indicesDists, searchParameters);
    Indices result(nFound);
    distanceOut.resize(nFound);
    for(UnsignedInteger k = 0; k < nFound; ++ k)
    {
      result[k] = indicesDists[k].first;
      distanceOut[k] = std::sqrt(indicesDists[k].second);
    }
    return result;
  }

private:
  Pointer<KDTreeSampleAdaptor> sampleAdaptor_;
  Pointer<nano_kd_tree_t> indexAdaptor_;
};

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
  if (!size)
    throw InvalidArgumentException(HERE) << "IntersectionMesher expected a non-empty collection";
  Mesh result(coll[0]);
  for (UnsignedInteger i = 1; i < coll.getSize(); ++ i)
    result = build2(result, coll[i]);
  return result;
}

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

  Sample vertices(0, dimension);
  CloudMesher cloudMesher;
  Collection<Indices> simplexColl;

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
      throw InternalException(HERE) << "dd_DDMatrix2Poly failed i1=" << i1;

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
        throw InternalException(HERE) << "dd_DDMatrix2Poly failed i2=" << i2;

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
    } // mesh2 simplices loop
    dd_FreeMatrix(h1);
    dd_FreePolyhedra(p1);

  } // mesh1 simplices loop

  // free cddlib objects
  dd_FreeMatrix(m1);
  dd_FreeMatrix(m2);
  dd_free_global_constants();

  // eliminate duplicate vertices
  const UnsignedInteger fullSize = vertices.getSize();
  Indices compressedVertexMap(fullSize);
  compressedVertexMap.fill();
  if (recompress_)
  {
    compressedVertexMap.fill(fullSize);
    Indices compressedIndices;
    const KDTree2 tree(vertices);
    const Scalar tolerance = 1e-12 * vertices.computeRange().norm();
    for (UnsignedInteger i = 0; i < fullSize; ++ i)
    {
      // check if already marked
      if (compressedVertexMap[i] < fullSize)
        continue;

      // retrieve the indices of unique points, beware the last few decimals can actually differ
      Point distance;
      const Indices nearest(tree.queryRadius(vertices[i], tolerance, distance));

      // mark the whole set of unique points (contains i index)
      const UnsignedInteger currentSize = compressedIndices.getSize();
      for (UnsignedInteger k = 0; k < nearest.getSize(); ++ k)
        compressedVertexMap[nearest[k]] = currentSize;

      // store the point
      compressedIndices.add(i);
    }
    LOGDEBUG(OSS() << "recompression fullSize=" << fullSize << " compressedSize=" << compressedIndices.getSize());
    vertices = vertices.select(compressedIndices);
  }

  // copy simplices
  IndicesCollection simplices(simplexColl.getSize(), dimension + 1);
  for (UnsignedInteger i = 0; i < simplexColl.getSize(); ++ i)
    for (UnsignedInteger j = 0; j <= dimension; ++ j)
      simplices(i, j) = compressedVertexMap[simplexColl[i][j]];

  return Mesh(vertices, simplices);
#else
  throw NotYetImplementedException(HERE) << "No cddlib support";
#endif
}


Mesh IntersectionMesher::buildConvex(const Collection<Mesh> & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (!size)
    throw InvalidArgumentException(HERE) << "IntersectionMesher expected a non-empty collection";
  Mesh result(coll[0]);
  for (UnsignedInteger i = 1; i < coll.getSize(); ++ i)
    result = build2Convex(result, coll[i]);
  return result;
}


Mesh IntersectionMesher::build2Convex(const Mesh & mesh1, const Mesh & mesh2) const
{
  const UnsignedInteger dimension = mesh1.getDimension();
  if (mesh2.getDimension() != dimension)
    throw InvalidArgumentException(HERE) << "IntersectionMesher expected meshes of same dimension";

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
