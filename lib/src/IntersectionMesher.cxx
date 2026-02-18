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
#include "otmeshing/ConvexDecompositionMesher.hxx"

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
      result = CompressMesh(result);
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

/* Deduplicate mesh vertices */
Mesh IntersectionMesher::CompressMesh(const Mesh & mesh)
{
  const UnsignedInteger dimension = mesh.getDimension();
  const Sample vertices(mesh.getVertices());
  const UnsignedInteger fullSize = vertices.getSize();
  if (!fullSize)
    return mesh;
  IndicesCollection simplices(mesh.getSimplices());
  Indices compressedVertexMap(fullSize, fullSize);
  const KDTree2 tree(vertices);
  const Scalar tolerance = SpecFunc::Precision * vertices.computeRange().norm();
  Sample verticesCompressed(0, dimension);
  for (UnsignedInteger i = 0; i < fullSize; ++ i)
  {
    // check if already marked
    if (compressedVertexMap[i] < fullSize)
      continue;

    // retrieve the indices of unique points, beware the last few decimals can actually differ
    Point distance;
    const Indices nearest(tree.queryRadius(vertices[i], tolerance, distance));

    // mark the whole set of unique points (contains i index)
    const UnsignedInteger currentSize = verticesCompressed.getSize();
    for (UnsignedInteger k = 0; k < nearest.getSize(); ++ k)
      compressedVertexMap[nearest[k]] = currentSize;

    // store the point
    verticesCompressed.add(vertices[i]);
  }
  LOGDEBUG(OSS() << "recompression fullSize=" << fullSize << " compressedSize=" << verticesCompressed.getSize());

  // renumber vertex indices
  for (UnsignedInteger i = 0; i < simplices.getSize(); ++ i)
    for (UnsignedInteger j = 0; j <= dimension; ++ j)
      simplices(i, j) = compressedVertexMap[simplices(i, j)];
  return Mesh(verticesCompressed, simplices);
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
          simplex.fill(verticesSize); // orientation may be incorrect
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

  Mesh result(vertices, IndicesCollection(simplexColl));
  if (recompress_)
    result = CompressMesh(result);
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
  Sample vertices(0, dimension);
  CloudMesher cloudMesher;
  Collection<Indices> simplexColl;
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
      simplex.fill(verticesSize); // orientation may be incorrect
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

  dd_FreeMatrix(gen);
  dd_free_global_constants();
  
  return Mesh(vertices, IndicesCollection(simplexColl));
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
