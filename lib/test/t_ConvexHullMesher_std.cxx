//                                               -*- C++ -*-
/**
 *  @brief The test file of class ConvexHullMesher
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
#include <cmath>

#include "otmeshing/otmeshing.hxx"
#include "otmeshing/ConvexHullMesher.hxx"
#include "openturns/OTtestcode.hxx"
#include "openturns/Box.hxx"
#include "openturns/Normal.hxx"

using namespace OT;
using namespace OTMESHING;
using namespace OT::Test;

int main()
{
  TESTPREAMBLE;

  ConvexHullMesher mesher;
  std::cout << mesher << std::endl;
  assert_equal(mesher.__repr__(), String("class=ConvexHullMesher"));

  // 1. empty sample → error
  try
  {
    mesher.build(Sample(0, 0));
    throw TestFailed(OSS() << "Expected InvalidArgumentException");
  }
  catch (const InvalidArgumentException &)
  {
    std::cout << "Empty sample correctly raised" << std::endl;
  }

  // 2. not enough points
  try
  {
    mesher.build(Sample(1, 1));
    throw TestFailed(OSS() << "Expected InvalidArgumentException");
  }
  catch (const InvalidArgumentException &)
  {
    std::cout << "Not enough points correctly raised" << std::endl;
  }

  // 3. 1D: segment [0.5, 1.5]
  {
    Sample points(0, 1);
    points.add(Point({0.5}));
    points.add(Point({1.5}));
    const Mesh hull(mesher.build(points));
    std::cout << "1D segment: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(1));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(1));
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(2));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(1));
    assert_almost_equal(hull.getVolume(), 1.0);
    assert_equal(hull.isValid(), true);
  }

  // 4. 1D: many collinear points → hull is [min, max]
  {
    Sample points(0, 1);
    points.add(Point({0.5}));
    points.add(Point({1.5}));
    points.add(Point({2.0}));
    points.add(Point({3.0}));
    points.add(Point({6.0}));
    const Mesh hull(mesher.build(points));
    std::cout << "1D collinear: " << hull << std::endl;
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(2));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(1));
    assert_almost_equal(hull.getVolume(), 5.5);
  }

  // 5. 2D: triangle
  {
    Sample points(0, 2);
    points.add(Point({0.0, 0.0}));
    points.add(Point({1.0, 0.0}));
    points.add(Point({0.0, 1.0}));
    const Mesh hull(mesher.build(points));
    std::cout << "2D triangle: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(2));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(1));
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(3));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(3));
    assert_almost_equal(hull.getVolume(), 2.0 + std::sqrt(2.0));
  }

  // 6. 2D: square
  {
    Sample points(0, 2);
    points.add(Point({0.0, 0.0}));
    points.add(Point({1.0, 0.0}));
    points.add(Point({1.0, 1.0}));
    points.add(Point({0.0, 1.0}));
    const Mesh hull(mesher.build(points));
    std::cout << "2D square: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(2));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(1));
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(hull.getVolume(), 4.0);
    assert_equal(hull.isValid(), true);
  }

  // 7. 2D: square with interior points → hull unchanged
  {
    Sample points(0, 2);
    points.add(Point({0.0, 0.0}));
    points.add(Point({1.0, 0.0}));
    points.add(Point({1.0, 1.0}));
    points.add(Point({0.0, 1.0}));
    points.add(Point({0.5, 0.5}));
    points.add(Point({0.2, 0.3}));
    const Mesh hull(mesher.build(points));
    std::cout << "2D square+interior: " << hull << std::endl;
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(hull.getVolume(), 4.0);
  }

  // 8. 3D: tetrahedron
  {
    Sample points(0, 3);
    points.add(Point({0.0, 0.0, 0.0}));
    points.add(Point({1.0, 0.0, 0.0}));
    points.add(Point({0.0, 1.0, 0.0}));
    points.add(Point({0.0, 0.0, 1.0}));
    const Mesh hull(mesher.build(points));
    std::cout << "3D tetrahedron: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(2));
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(4));
    // 4 triangular faces: areas are 0.5, 0.5, 0.5, sqrt(3)/2
    assert_almost_equal(hull.getVolume(), 1.5 + 0.5 * std::sqrt(3.0));
    assert_equal(hull.isValid(), true);
  }

  // 9. 3D: cube — critical test (triggers "Qt" centroid in Qhull)
  {
    Sample points(0, 3);
    points.add(Point({0.0, 0.0, 0.0}));
    points.add(Point({1.0, 0.0, 0.0}));
    points.add(Point({1.0, 1.0, 0.0}));
    points.add(Point({0.0, 1.0, 0.0}));
    points.add(Point({0.0, 0.0, 1.0}));
    points.add(Point({1.0, 0.0, 1.0}));
    points.add(Point({1.0, 1.0, 1.0}));
    points.add(Point({0.0, 1.0, 1.0}));
    const Mesh hull(mesher.build(points));
    std::cout << "3D cube: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(2));
    // CGAL: 8 vertices; Qhull: 14 (8 + 6 face centroids)
    if (hull.getVerticesNumber() < 8)
      throw TestFailed(OSS() << "Expected at least 8 vertices, got " << hull.getVerticesNumber());
    // 6 faces × 2 triangles = 12
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(12));
    assert_almost_equal(hull.getVolume(), 6.0);
    assert_equal(hull.isValid(), true);
  }

  // 10. 4D: simplex (5 vertices)
  {
    Sample points(0, 4);
    points.add(Point({0.0, 0.0, 0.0, 0.0}));
    points.add(Point({1.0, 0.0, 0.0, 0.0}));
    points.add(Point({0.0, 1.0, 0.0, 0.0}));
    points.add(Point({0.0, 0.0, 1.0, 0.0}));
    points.add(Point({0.0, 0.0, 0.0, 1.0}));
    const Mesh hull(mesher.build(points));
    std::cout << "4D simplex: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(4));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(hull.getVerticesNumber(), static_cast<UnsignedInteger>(5));
    assert_equal(hull.getSimplicesNumber(), static_cast<UnsignedInteger>(5));
    assert_almost_equal(hull.getVolume(), 1.0);
    assert_equal(hull.isValid(), true);
  }

  // 11. 4D: hypercube (16 vertices)
  {
    const Box box(Indices(4, 0));
    const Sample points(box.generate());
    const Mesh hull(mesher.build(points));
    std::cout << "4D cube: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(4));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(3));
    if (hull.getVerticesNumber() < 16)
      throw TestFailed(OSS() << "Expected at least 16 vertices, got " << hull.getVerticesNumber());
    assert_almost_equal(hull.getVolume(), 8.0);
    assert_equal(hull.isValid(), true);
  }

  // 12. reusing the same mesher for multiple builds
  {
    Sample points1(0, 2);
    points1.add(Point({0.0, 0.0}));
    points1.add(Point({1.0, 0.0}));
    points1.add(Point({0.0, 1.0}));
    const Mesh hull1(mesher.build(points1));
    assert_equal(hull1.getVerticesNumber(), static_cast<UnsignedInteger>(3));

    Sample points2(0, 2);
    points2.add(Point({0.0, 0.0}));
    points2.add(Point({2.0, 0.0}));
    points2.add(Point({0.0, 2.0}));
    const Mesh hull2(mesher.build(points2));
    assert_equal(hull2.getVerticesNumber(), static_cast<UnsignedInteger>(3));
  }

  // 13. Gaussian random sample (dim 2, spot check)
  {
    const Sample points(Normal(2).getSample(500));
    const Mesh hull(mesher.build(points));
    std::cout << "2D gaussian: " << hull << std::endl;
    assert_equal(hull.getDimension(), static_cast<UnsignedInteger>(2));
    assert_equal(hull.getIntrinsicDimension(), static_cast<UnsignedInteger>(1));
    if (hull.getVerticesNumber() >= points.getSize())
      throw TestFailed(OSS() << "Expected fewer hull vertices than input points");
    assert_equal(hull.isValid(), true);
  }

  return 0;
}
