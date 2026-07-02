//                                               -*- C++ -*-
/**
 *  @brief The test file of class IntersectionMesher
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
#include <cmath>
#include <cstdio>
#include <iostream>

#include <openturns/OT.hxx>
#include <openturns/OTtestcode.hxx>

#include "otmeshing/otmeshing.hxx"

using namespace OT;
using namespace OTMESHING;
using namespace OT::Test;

typedef IntersectionMesher::MeshCollection MeshCollection;
typedef IntersectionMesher::SampleCollection SampleCollection;
typedef IntersectionMesher::CylinderCollection CylinderCollection;

int main()
{
  TESTPREAMBLE;

  // 1. Default construction and string representation
  IntersectionMesher mesher;
  std::cout << mesher << std::endl;
  std::cout << "recompress=" << mesher.getRecompress() << std::endl;
  std::cout << "useSimplicesDecomposition=" << mesher.getUseSimplicesDecomposition() << std::endl;

  // 2. Empty collection -> empty mesh
  {
    const Mesh empty(mesher.build(MeshCollection()));
    std::cout << "Empty collection: " << empty << std::endl;
    assert_equal(empty.getDimension(), 0UL);
    assert_equal(empty.getVerticesNumber(), 0UL);
  }

  // 3. Single mesh -> returns it
  for (UnsignedInteger dim = 2; dim < 5; ++dim)
  {
    const Mesh mesh(IntervalMesher(Indices(dim, 1)).build(Interval(dim)));
    const Mesh single(mesher.build(MeshCollection(1, mesh)));
    std::cout << "Single dim=" << dim << ": " << single << std::endl;
    assert_equal(single.getVerticesNumber(), mesh.getVerticesNumber());
    assert_equal(single.getSimplicesNumber(), mesh.getSimplicesNumber());
    assert_almost_equal(single.getVolume(), 1.0);
  }

  // 4. Overlapping cubes (2D-5D)
  for (UnsignedInteger dim = 2; dim < 6; ++dim)
  {
    const Mesh mesh1(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 0.0), Point(dim, 3.0))));
    const Mesh mesh2(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 1.0), Point(dim, 4.0))));
    const Mesh inter(mesher.build(MeshCollection({mesh1, mesh2})));
    std::cout << "Overlap dim=" << dim << ": " << inter << " vol=" << inter.getVolume() << std::endl;
    assert_almost_equal(inter.getVolume(), std::pow(2.0, (int)dim));
  }

  // 5. Disjoint cubes -> empty intersection
  for (UnsignedInteger dim = 2; dim < 5; ++dim)
  {
    const Mesh mesh1(IntervalMesher(Indices(dim, 1)).build(Interval(dim)));
    const Mesh mesh2(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 3.0), Point(dim, 4.0))));
    const Mesh inter(mesher.build(MeshCollection({mesh1, mesh2})));
    std::cout << "Disjoint dim=" << dim << ": " << inter << std::endl;
    assert_equal(inter.getVerticesNumber(), 0UL);
  }

  // 6. Self-intersection (mesh with itself)
  {
    const Mesh mesh(IntervalMesher(Indices(3, 2)).build(Interval(Point(3, 0.0), Point(3, 3.0))));
    const Mesh inter(mesher.build(MeshCollection({mesh, mesh})));
    std::cout << "Self-intersection: " << inter << std::endl;
    assert_almost_equal(inter.getVolume(), mesh.getVolume());
  }

  // 7. Three-mesh intersection (nested)
  {
    const Mesh mesh1(IntervalMesher(Indices(3, 2)).build(Interval(Point(3, 0.0), Point(3, 3.0))));
    const Mesh mesh2(IntervalMesher(Indices(3, 2)).build(Interval(Point(3, 0.5), Point(3, 2.5))));
    const Mesh mesh3(IntervalMesher(Indices(3, 2)).build(Interval(Point(3, 1.0), Point(3, 2.0))));
    const Mesh inter(mesher.build(MeshCollection({mesh1, mesh2, mesh3})));
    std::cout << "Three-mesh nested: " << inter << std::endl;
    assert_almost_equal(inter.getVolume(), 1.0);
  }

  // 8. buildConvex with overlapping convex meshes (2D-5D)
  for (UnsignedInteger dim = 2; dim < 6; ++dim)
  {
    const Mesh mesh1(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 0.0), Point(dim, 3.0))));
    const Mesh mesh2(IntervalMesher(Indices(dim, 1)).build(Interval(Point(dim, 1.0), Point(dim, 4.0))));
    const Mesh inter(mesher.buildConvex(MeshCollection({mesh1, mesh2})));
    std::cout << "buildConvex dim=" << dim << ": " << inter << " vol=" << inter.getVolume() << std::endl;
    assert_almost_equal(inter.getVolume(), std::pow(2.0, (int)dim));
  }

  // 9. buildConvexSample with two samples
  {
    SampleCollection coll;
    {
      Sample s1(0, 2);
      s1.add(Point({0.0, 0.0}));
      s1.add(Point({1.0, 0.0}));
      s1.add(Point({1.0, 1.0}));
      s1.add(Point({0.0, 1.0}));
      coll.add(s1);
    }
    {
      Sample s2(0, 2);
      s2.add(Point({0.5, 0.5}));
      s2.add(Point({1.5, 0.5}));
      s2.add(Point({1.5, 1.5}));
      s2.add(Point({0.5, 1.5}));
      coll.add(s2);
    }
    const Sample result(mesher.buildConvexSample(coll));
    std::cout << "buildConvexSample pairwise: " << result << std::endl;
    OT::Test::assert(result.getSize() >= 3);
  }

  // 10. buildConvexSample two-argument overload
  {
    Sample s1(0, 2);
    s1.add(Point({0.0, 0.0}));
    s1.add(Point({1.0, 0.0}));
    s1.add(Point({1.0, 1.0}));
    s1.add(Point({0.0, 1.0}));
    Sample s2(0, 2);
    s2.add(Point({0.5, 0.5}));
    s2.add(Point({1.5, 0.5}));
    s2.add(Point({1.5, 1.5}));
    s2.add(Point({0.5, 1.5}));
    const Sample result(mesher.buildConvexSample(s1, s2));
    std::cout << "buildConvexSample two-arg: " << result << std::endl;
    OT::Test::assert(result.getSize() >= 3);
  }

  // 11. buildConvexSample dimension mismatch should throw
  {
    Sample s1(0, 2);
    s1.add(Point({0.0, 0.0}));
    s1.add(Point({1.0, 0.0}));
    Sample s2(0, 3);
    s2.add(Point({0.0, 0.0, 0.0}));
    s2.add(Point({1.0, 0.0, 0.0}));
    try
    {
      mesher.buildConvexSample(s1, s2);
      throw TestFailed(OSS() << "Expected InvalidArgumentException");
    }
    catch (const InvalidArgumentException &)
    {
      std::cout << "Dimension mismatch correctly raised for buildConvexSample" << std::endl;
    }
  }

  // 12. buildConvex empty / single / self
  {
    const Mesh empty(mesher.buildConvex(MeshCollection()));
    assert_equal(empty.getDimension(), 0UL);

    const Mesh mesh(IntervalMesher(Indices(3, 1)).build(Interval(3)));
    const Mesh single(mesher.buildConvex(MeshCollection({mesh})));
    assert_almost_equal(single.getVolume(), mesh.getVolume());

    const Mesh self(mesher.buildConvex(MeshCollection({mesh, mesh})));
    assert_almost_equal(self.getVolume(), mesh.getVolume());
  }

  // 13. buildConvexSample empty / single collections
  {
    const Sample empty(mesher.buildConvexSample(SampleCollection()));
    assert_equal(empty.getSize(), 0UL);

    Sample s1(0, 2);
    s1.add(Point({0.0, 0.0}));
    s1.add(Point({1.0, 0.0}));
    s1.add(Point({1.0, 1.0}));
    s1.add(Point({0.0, 1.0}));
    const Sample single(mesher.buildConvexSample(SampleCollection({s1})));
    assert_equal(single.getSize(), s1.getSize());
  }

  // 14. buildCylinder: convex cylinder intersection (Steinmetz solid)
  {
    const UnsignedInteger nTheta = 32;
    const Scalar R = 2.0;
    const Scalar H = 5.0;
    Sample circle1(0, 2);
    Sample circle2(0, 2);
    for (UnsignedInteger i = 0; i < nTheta; ++i)
    {
      const Scalar theta = i * 2.0 * M_PI / nTheta;
      circle1.add(Point({R * std::cos(theta), R * std::sin(theta)}));
      circle2.add(Point({R * std::cos(theta), R * std::sin(theta)}));
    }
    const Mesh disc1(PolygonMesher().build(circle1));
    const Mesh disc2(PolygonMesher().build(circle2));
    const Cylinder cyl1(disc1, Interval(Point({-H/2}), Point({H/2})), Indices({2}), 2);
    const Cylinder cyl2(disc2, Interval(Point({-H/2}), Point({H/2})), Indices({0}), 2);
    const Mesh inter(mesher.buildCylinder(CylinderCollection({cyl1, cyl2})));
    std::cout << "Cylinder intersection vol=" << inter.getVolume() << std::endl;
    assert_almost_equal(inter.getVolume(), 16.0 / 3.0 * R * R * R, 2e-1);
  }

  // 15. buildCylinder: non-convex cylinder intersection
  {
    const UnsignedInteger nTheta = 32;
    const Scalar R = 2.0;
    const Scalar H = 5.0;
    Sample circle(0, 2);
    Sample star(0, 2);
    for (UnsignedInteger i = 0; i < nTheta; ++i)
    {
      const Scalar theta = i * 2.0 * M_PI / nTheta;
      const Scalar x1 = R * std::cos(theta);
      const Scalar y1 = R * std::sin(theta);
      circle.add(Point({x1, y1}));
      if (i % 2 == 0)
        star.add(Point({x1, y1}));
      else
        star.add(Point({x1 * 1.5, y1 * 1.5}));
    }
    const Mesh disc2(PolygonMesher().build(circle));
    const Mesh disc3(PolygonMesher().build(star));
    const Cylinder cyl2(disc2, Interval(Point({-H/2}), Point({H/2})), Indices({0}), 2);
    const Cylinder cyl3(disc3, Interval(Point({-H/2}), Point({H/2})), Indices({2}), 2);
    const Mesh inter(mesher.buildCylinder(CylinderCollection({cyl3, cyl2})));
    std::cout << "Non-convex cylinder intersection vol=" << inter.getVolume() << std::endl;
    assert_almost_equal(inter.getVolume(), 53.4976, 1e-2);
  }

  // 16. buildCylinder: disjoint cylinders -> empty
  {
    const Mesh disc(IntervalMesher(Indices(2, 2)).build(Interval(Point({-20.0, -20.0}), Point({-10.0, -10.0}))));
    const Cylinder cyl(disc, Interval(Point({-50.0}), Point({-40.0})), Indices({2}), 2);
    const Cylinder cyl2(disc, Interval(Point({-50.0}), Point({-40.0})), Indices({0}), 2);
    const Mesh inter(mesher.buildCylinder(CylinderCollection({cyl, cyl2})));
    std::cout << "Disjoint cylinders: " << inter << std::endl;
    assert_equal(inter.getVerticesNumber(), 0UL);
  }

  // 17. buildCylinder: self-intersection
  {
    const Mesh disc(IntervalMesher(Indices(2, 2)).build(Interval(Point({-1.0, -1.0}), Point({1.0, 1.0}))));
    const Cylinder cyl(disc, Interval(Point({-1.0}), Point({1.0})), Indices({2}), 2);
    const Mesh inter(mesher.buildCylinder(CylinderCollection({cyl, cyl})));
    std::cout << "Cylinder self-intersection vol=" << inter.getVolume() << std::endl;
    const Mesh expected(cyl.computeMesh());
    assert_almost_equal(inter.getVolume(), expected.getVolume());
  }

  // 18. Recompress flag
  {
    IntersectionMesher m;
    assert_equal(m.getRecompress(), true);
    m.setRecompress(false);
    assert_equal(m.getRecompress(), false);
    m.setRecompress(true);
    assert_equal(m.getRecompress(), true);
  }

  // 19. UseSimplicesDecomposition flag
  {
    IntersectionMesher m;
    assert_equal(m.getUseSimplicesDecomposition(), true);
    m.setUseSimplicesDecomposition(false);
    assert_equal(m.getUseSimplicesDecomposition(), false);
    m.setUseSimplicesDecomposition(true);
    assert_equal(m.getUseSimplicesDecomposition(), true);
  }

  // 20. Save/Load
  {
    const String studyFile("t_IntersectionMesher_std.xml");
    IntersectionMesher m1;
    m1.setRecompress(false);
    m1.setUseSimplicesDecomposition(false);
    {
      Study study(studyFile);
      study.add("mesher", m1);
      study.save();
    }
    {
      IntersectionMesher m2;
      Study loadedStudy(studyFile);
      loadedStudy.load();
      loadedStudy.fillObject("mesher", m2);
      std::cout << "Loaded: recompress=" << m2.getRecompress()
                << " useSimplices=" << m2.getUseSimplicesDecomposition() << std::endl;
      assert_equal(m2.getRecompress(), false);
      assert_equal(m2.getUseSimplicesDecomposition(), false);
    }
    std::remove(studyFile.c_str());
  }

  return 0;
}
