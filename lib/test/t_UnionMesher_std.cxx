//                                               -*- C++ -*-
/**
 *  @brief The test file of class UnionMesher
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
#include <iostream>

#include "otmeshing/otmeshing.hxx"
#include "openturns/OTtestcode.hxx"
#include "openturns/IntervalMesher.hxx"

using namespace OT;
using namespace OTMESHING;
using namespace OT::Test;

typedef UnionMesher::MeshCollection MeshCollection;

int main()
{
  TESTPREAMBLE;

  // 1. Default construction and string representation
  UnionMesher mesher;
  std::cout << mesher << std::endl;

  // 2. Empty collection → empty mesh (dimension 0)
  {
    Mesh empty(mesher.build(MeshCollection()));
    std::cout << "Empty collection: " << empty << std::endl;
    assert_equal(empty.getDimension(), static_cast<UnsignedInteger>(0));
    assert_equal(empty.getVerticesNumber(), static_cast<UnsignedInteger>(0));
  }

  // 3. Single mesh → returns compressed mesh (no-op for valid input)
  for (UnsignedInteger dim = 2; dim < 5; ++dim)
  {
    const Mesh mesh(IntervalMesher(Indices(dim, 1)).build(Interval(dim)));
    const Mesh single(mesher.build(MeshCollection(1, mesh)));
    std::cout << "Single dim=" << dim << ": " << single << std::endl;
    assert_equal(single.getVerticesNumber(), mesh.getVerticesNumber());
    assert_equal(single.getSimplicesNumber(), mesh.getSimplicesNumber());
    assert_almost_equal(single.getVolume(), 1.0);
  }

  // 4. Non-overlapping meshes
  for (UnsignedInteger dim = 2; dim < 6; ++dim)
  {
    const Mesh mesh1(IntervalMesher(Indices(dim, 1)).build(Interval(dim)));
    const Mesh mesh2(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 2.0), Point(dim, 3.0))));
    const Mesh unionMesh(mesher.build(MeshCollection({mesh1, mesh2})));
    std::cout << "Disjoint dim=" << dim << ": " << unionMesh << std::endl;
    assert_almost_equal(unionMesh.getVolume(), 2.0);
  }

  // 5. Touching meshes (shared boundary) — exercises vertex dedup
  // mesh1 = [0,1]^dim, mesh2 = [1,2]x[0,1]^{dim-1}, sharing the face at x=1
  for (UnsignedInteger dim = 2; dim < 5; ++dim)
  {
    const Mesh mesh1(IntervalMesher(Indices(dim, 1)).build(Interval(dim)));
    Point lower(dim, 0.0);
    Point upper(dim, 1.0);
    lower[0] = 1.0;
    upper[0] = 2.0;
    const Mesh mesh2(IntervalMesher(Indices(dim, 1)).build(Interval(lower, upper)));
    const Mesh unionMesh(mesher.build(MeshCollection({mesh1, mesh2})));
    std::cout << "Touching dim=" << dim << ": " << unionMesh << std::endl;
    assert_almost_equal(unionMesh.getVolume(), 2.0);
    // shared (dim-1)-face has 2^{dim-1} vertices
    const UnsignedInteger nNoDedup = (1u << (dim + 1));
    const UnsignedInteger nShared = (1u << (dim - 1));
    assert_equal(unionMesh.getVerticesNumber(), nNoDedup - nShared);
  }

  // 6. Three disjoint meshes
  for (UnsignedInteger dim = 2; dim < 4; ++dim)
  {
    const Mesh mesh1(IntervalMesher(Indices(dim, 1)).build(Interval(dim)));
    const Mesh mesh2(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 2.0), Point(dim, 3.0))));
    const Mesh mesh3(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 4.0), Point(dim, 5.0))));
    const Mesh unionMesh(mesher.build(MeshCollection({mesh1, mesh2, mesh3})));
    std::cout << "Three dim=" << dim << ": " << unionMesh << std::endl;
    assert_almost_equal(unionMesh.getVolume(), 3.0);
  }

  // 7. CompressMesh static method: no-op on clean mesh
  {
    const Mesh mesh(IntervalMesher(Indices(2, 1)).build(Interval(2)));
    const Mesh compressed(UnionMesher::CompressMesh(mesh));
    std::cout << "Already compressed: " << compressed << std::endl;
    assert_equal(compressed.getVerticesNumber(), mesh.getVerticesNumber());
    assert_equal(compressed.getSimplicesNumber(), mesh.getSimplicesNumber());
    assert_almost_equal(compressed.getVolume(), 1.0);
  }

  // 8. CompressMesh: collapse duplicate vertices
  {
    Sample vertices(0, 2);
    vertices.add(Point({0.0, 0.0}));
    vertices.add(Point({1.0, 0.0}));
    vertices.add(Point({1.0, 1.0}));
    vertices.add(Point({0.0, 1.0}));
    vertices.add(Point({0.0, 0.0}));
    vertices.add(Point({1.0, 0.0}));
    IndicesCollection simplices(4, 3);
    simplices(0, 0) = 0; simplices(0, 1) = 1; simplices(0, 2) = 2;
    simplices(1, 0) = 0; simplices(1, 1) = 2; simplices(1, 2) = 3;
    simplices(2, 0) = 4; simplices(2, 1) = 5; simplices(2, 2) = 2;
    simplices(3, 0) = 4; simplices(3, 1) = 2; simplices(3, 2) = 3;
    Mesh mesh(vertices, simplices, false);
    const Mesh compressed(UnionMesher::CompressMesh(mesh));
    std::cout << "Duplicates collapsed: " << compressed << std::endl;
    assert_equal(compressed.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(compressed.getVolume(), 2.0);
  }

  // 9. CompressMesh: unused vertices dropped
  {
    Sample vertices(0, 2);
    vertices.add(Point({0.0, 0.0}));
    vertices.add(Point({1.0, 0.0}));
    vertices.add(Point({1.0, 1.0}));
    vertices.add(Point({0.0, 1.0}));
    vertices.add(Point({42.0, 42.0}));
    IndicesCollection simplices(2, 3);
    simplices(0, 0) = 0; simplices(0, 1) = 1; simplices(0, 2) = 2;
    simplices(1, 0) = 0; simplices(1, 1) = 2; simplices(1, 2) = 3;
    Mesh mesh(vertices, simplices, false);
    const Mesh compressed(UnionMesher::CompressMesh(mesh));
    std::cout << "Unused dropped: " << compressed << std::endl;
    assert_equal(compressed.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(compressed.getVolume(), 1.0);
  }

  // 10. build with a mesh that has an unused vertex
  {
    Sample vertices(0, 2);
    vertices.add(Point({0.0, 0.0}));
    vertices.add(Point({1.0, 0.0}));
    vertices.add(Point({1.0, 1.0}));
    vertices.add(Point({0.0, 1.0}));
    vertices.add(Point({42.0, 42.0}));
    IndicesCollection simplices(2, 3);
    simplices(0, 0) = 0; simplices(0, 1) = 1; simplices(0, 2) = 2;
    simplices(1, 0) = 0; simplices(1, 1) = 2; simplices(1, 2) = 3;
    Mesh mesh(vertices, simplices, false);
    const Mesh unionMesh(mesher.build(MeshCollection({mesh})));
    std::cout << "Build with unused: " << unionMesh << std::endl;
    assert_equal(unionMesh.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(unionMesh.getVolume(), 1.0);
  }

  // 11. Empty input mesh
  {
    const Mesh emptyMesh(Sample(0, 2), IndicesCollection());
    const Mesh single(mesher.build(MeshCollection({emptyMesh})));
    std::cout << "Empty mesh: " << single << std::endl;
  }

  // 12. Two empty meshes
  {
    const Mesh emptyMesh(Sample(0, 2), IndicesCollection());
    const Mesh unionMesh(mesher.build(MeshCollection({emptyMesh, emptyMesh})));
    std::cout << "Two empty: " << unionMesh << std::endl;
  }

  // 13. Dimension mismatch should throw
  {
    const Mesh mesh1(IntervalMesher(Indices(2, 1)).build(Interval(2)));
    const Mesh mesh2(IntervalMesher(Indices(3, 1)).build(Interval(3)));
    try
    {
      mesher.build(MeshCollection({mesh1, mesh2}));
      throw TestFailed(OSS() << "Expected InvalidArgumentException");
    }
    catch (const InvalidArgumentException &)
    {
      std::cout << "Dimension mismatch correctly raised" << std::endl;
    }
  }

  // 14. CompressMesh: duplicate vertices with non-zero range (tolerance > 0)
  {
    Sample vertices(0, 2);
    vertices.add(Point({0.0, 0.0}));
    vertices.add(Point({1.0, 0.0}));
    vertices.add(Point({1.0, 1.0}));
    vertices.add(Point({0.0, 1.0}));
    vertices.add(Point({0.0, 0.0}));
    vertices.add(Point({1.0, 0.0}));
    IndicesCollection simplices(4, 3);
    simplices(0, 0) = 0; simplices(0, 1) = 1; simplices(0, 2) = 2;
    simplices(1, 0) = 0; simplices(1, 1) = 2; simplices(1, 2) = 3;
    simplices(2, 0) = 4; simplices(2, 1) = 5; simplices(2, 2) = 2;
    simplices(3, 0) = 4; simplices(3, 1) = 2; simplices(3, 2) = 3;
    Mesh mesh(vertices, simplices, false);
    const Mesh compressed(UnionMesher::CompressMesh(mesh));
    std::cout << "Duplicates non-zero range: " << compressed << std::endl;
    assert_equal(compressed.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(compressed.getVolume(), 2.0);
  }

  return 0;
}
