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

Mesh UnionMesher::build(const MeshCollection & coll) const
{
  const UnsignedInteger size = coll.getSize();
  if (size == 0)
    return Mesh(Sample(0, 0));
  else if (size == 1)
    return coll[0];

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
  return Mesh(vertices, simplices);
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
