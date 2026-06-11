//                                               -*- C++ -*-
/**
 *  @brief The test file of class ConvexDecompositionMesher
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

#include <openturns/OT.hxx>
#include <openturns/OTtestcode.hxx>

#include "otmeshing/ConvexDecompositionMesher.hxx"
#include "otmeshing/PolygonMesher.hxx"
#include "otmeshing/UnionMesher.hxx"

using namespace OT;
using namespace OT::Test;
using namespace OTMESHING;

int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {
    // Default constructor
    fullprint << "--- Default constructor ---" << std::endl;
    ConvexDecompositionMesher mesher;
    assert_equal(mesher.getUseSimplicesDecomposition(), false);
    fullprint << "OK" << std::endl;

    // __repr__
    fullprint << "--- __repr__ ---" << std::endl;
    assert_equal(mesher.__repr__(), String("class=ConvexDecompositionMesher"));
    fullprint << "OK" << std::endl;

    // Flag roundtrip
    fullprint << "--- useSimplicesDecomposition get/set ---" << std::endl;
    mesher.setUseSimplicesDecomposition(true);
    assert_equal(mesher.getUseSimplicesDecomposition(), true);
    mesher.setUseSimplicesDecomposition(false);
    assert_equal(mesher.getUseSimplicesDecomposition(), false);
    fullprint << "OK" << std::endl;

    // IsConvex on convex 2D mesh
    fullprint << "--- IsConvex (convex 2D mesh) ---" << std::endl;
    const Mesh convexMesh2D(IntervalMesher(Indices({2, 2})).build(Interval(Point({0.0, 0.0}), Point({1.0, 1.0}))));
    assert_equal(ConvexDecompositionMesher::IsConvex(convexMesh2D), true);
    fullprint << "OK" << std::endl;

    // IsConvex on convex 3D mesh
    fullprint << "--- IsConvex (convex 3D mesh) ---" << std::endl;
    const Mesh convexMesh3D(IntervalMesher(Indices({1, 1, 1})).build(Interval(Point({0.0, 0.0, 0.0}), Point({1.0, 1.0, 1.0}))));
    assert_equal(ConvexDecompositionMesher::IsConvex(convexMesh3D), true);
    fullprint << "OK" << std::endl;

    // 2D build on convex mesh — volume sum invariant
    fullprint << "--- build 2D convex mesh (volume invariant) ---" << std::endl;
    {
      Collection<Mesh> decomp = mesher.build(convexMesh2D);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      assert_almost_equal(volumeSum, convexMesh2D.getVolume());
      fullprint << "Number of components: " << decomp.getSize() << ", volume: " << volumeSum << ", OK" << std::endl;
    }

    // 3D build on convex mesh returns 1 component
    fullprint << "--- build 3D convex -> 1 component ---" << std::endl;
    {
      Collection<Mesh> decomp = mesher.build(convexMesh3D);
      assert_equal(decomp.getSize(), UnsignedInteger(1));
      assert_equal(ConvexDecompositionMesher::IsConvex(decomp[0]), true);
      fullprint << "OK" << std::endl;
    }

    // 2D build on non-convex mesh
    fullprint << "--- build 2D non-convex snake+box ---" << std::endl;
    {
      Sample polyline(0, 2);
      polyline.add(Point({0.0, 0.0}));
      polyline.add(Point({0.0, 5.0}));
      polyline.add(Point({6.0, 5.0}));
      polyline.add(Point({6.0, 0.0}));
      polyline.add(Point({2.0, 0.0}));
      polyline.add(Point({2.0, 3.0}));
      polyline.add(Point({4.0, 3.0}));
      polyline.add(Point({4.0, 2.0}));
      polyline.add(Point({3.0, 2.0}));
      polyline.add(Point({3.0, 1.0}));
      polyline.add(Point({5.0, 1.0}));
      polyline.add(Point({5.0, 4.0}));
      polyline.add(Point({1.0, 4.0}));
      polyline.add(Point({1.0, 0.0}));
      const Mesh snakeMesh = PolygonMesher().build(polyline);
      const Mesh boxMesh = IntervalMesher(Indices({1, 1})).build(Interval(Point({-2.0, -2.0}), Point({-1.0, -1.0})));
      const Mesh mesh = UnionMesher().build(Collection<Mesh>({snakeMesh, boxMesh}));
      assert_equal(ConvexDecompositionMesher::IsConvex(mesh), false);

      Collection<Mesh> decomp = mesher.build(mesh);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      assert_almost_equal(volumeSum, mesh.getVolume());
      fullprint << "Number of components: " << decomp.getSize() << ", OK" << std::endl;
    }

    // 3D build on non-convex surface mesh (notched cube)
    fullprint << "--- build 3D non-convex surface mesh ---" << std::endl;
    {
      Sample vertices(0, 3);
      vertices.add(Point({0.0, 0.0, 0.0}));
      vertices.add(Point({2.0, 0.0, 0.0}));
      vertices.add(Point({2.0, 2.0, 0.0}));
      vertices.add(Point({0.0, 2.0, 0.0}));
      vertices.add(Point({0.0, 0.0, 2.0}));
      vertices.add(Point({2.0, 0.0, 2.0}));
      vertices.add(Point({0.0, 2.0, 2.0}));
      vertices.add(Point({1.0, 1.0, 1.0}));
      vertices.add(Point({2.0, 1.0, 1.0}));
      vertices.add(Point({2.0, 2.0, 1.0}));
      vertices.add(Point({1.0, 2.0, 1.0}));
      vertices.add(Point({1.0, 1.0, 2.0}));
      vertices.add(Point({2.0, 1.0, 2.0}));
      vertices.add(Point({1.0, 2.0, 2.0}));

      Collection<Indices> simplicesColl;
      simplicesColl.add(Indices({0, 1, 5, 5}));
      simplicesColl.add(Indices({0, 5, 4, 4}));
      simplicesColl.add(Indices({0, 4, 6, 6}));
      simplicesColl.add(Indices({0, 6, 3, 3}));
      simplicesColl.add(Indices({0, 2, 1, 1}));
      simplicesColl.add(Indices({0, 3, 2, 2}));
      simplicesColl.add(Indices({10, 9, 2, 2}));
      simplicesColl.add(Indices({10, 2, 3, 3}));
      simplicesColl.add(Indices({10, 3, 6, 6}));
      simplicesColl.add(Indices({10, 6, 13, 13}));
      simplicesColl.add(Indices({8, 2, 9, 9}));
      simplicesColl.add(Indices({8, 1, 2, 2}));
      simplicesColl.add(Indices({8, 5, 1, 1}));
      simplicesColl.add(Indices({8, 12, 5, 5}));
      simplicesColl.add(Indices({11, 13, 6, 6}));
      simplicesColl.add(Indices({11, 6, 4, 4}));
      simplicesColl.add(Indices({11, 4, 5, 5}));
      simplicesColl.add(Indices({11, 5, 12, 12}));
      simplicesColl.add(Indices({8, 7, 12, 12}));
      simplicesColl.add(Indices({7, 11, 12, 12}));
      simplicesColl.add(Indices({7, 10, 13, 13}));
      simplicesColl.add(Indices({7, 13, 11, 11}));
      simplicesColl.add(Indices({7, 8, 9, 9}));
      simplicesColl.add(Indices({7, 9, 10, 10}));
      const IndicesCollection simplices(simplicesColl);

      const Mesh mesh(vertices, simplices);
      assert(mesh.isValid());
      assert_equal(ConvexDecompositionMesher::IsConvex(mesh), false);

      Collection<Mesh> decomp = mesher.build(mesh);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      // The notched box has volume = outer box (2*2*2=8) - inner notch (1*1*1=1) = 7
      assert_almost_equal(volumeSum, 7.0);
      fullprint << "Number of components: " << decomp.getSize() << ", volume: " << volumeSum << ", OK" << std::endl;
    }

    // 3D volumetric (two overlapping cubes)
    fullprint << "--- build 3D volumetric overlapping cubes ---" << std::endl;
    {
      const Mesh cube1Mesh = IntervalMesher(Indices({1, 1, 1})).build(Interval(Point({0.0, 0.0, 0.0}), Point({2.0, 2.0, 2.0})));
      const Mesh cube2Mesh = IntervalMesher(Indices({1, 1, 1})).build(Interval(Point({1.0, 1.0, 1.0}), Point({3.0, 3.0, 3.0})));
      const Mesh mesh = UnionMesher().build(Collection<Mesh>({cube1Mesh, cube2Mesh}));
      assert(mesh.isValid());

      Collection<Mesh> decomp = mesher.build(mesh);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      // Overlap: each cube is 2^3=8, total 16, overlapping region 1^3=1, so union = 15
      assert_almost_equal(volumeSum, 15.0);
      fullprint << "Number of components: " << decomp.getSize() << ", volume: " << volumeSum << ", OK" << std::endl;
    }

    // Simplices decomposition mode in 2D
    fullprint << "--- build 2D simplices decomposition ---" << std::endl;
    {
      ConvexDecompositionMesher simpleMesher;
      simpleMesher.setUseSimplicesDecomposition(true);

      Sample polyline(0, 2);
      polyline.add(Point({0.0, 0.0}));
      polyline.add(Point({2.0, 0.0}));
      polyline.add(Point({2.0, 1.0}));
      polyline.add(Point({0.0, 1.0}));
      const Mesh mesh = PolygonMesher().build(polyline);
      assert(mesh.isValid());

      Collection<Mesh> decomp = simpleMesher.build(mesh);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        // Each component should be a single simplex (triangle)
        assert_equal(decomp[i].getSimplicesNumber(), UnsignedInteger(1));
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      assert_almost_equal(volumeSum, mesh.getVolume());
      fullprint << "Number of triangles: " << decomp.getSize() << ", volume: " << volumeSum << ", OK" << std::endl;
    }

    // Simplices decomposition mode in 3D
    fullprint << "--- build 3D volumetric simplices decomposition ---" << std::endl;
    {
      ConvexDecompositionMesher simpleMesher;
      simpleMesher.setUseSimplicesDecomposition(true);

      const Mesh mesh = IntervalMesher(Indices({1, 1, 1})).build(Interval(Point({0.0, 0.0, 0.0}), Point({1.0, 1.0, 1.0})));
      assert(mesh.isValid());

      Collection<Mesh> decomp = simpleMesher.build(mesh);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      assert_almost_equal(volumeSum, mesh.getVolume());
      fullprint << "Number of components: " << decomp.getSize() << ", volume: " << volumeSum << ", OK" << std::endl;
    }

    // Clone
    fullprint << "--- clone ---" << std::endl;
    {
      ConvexDecompositionMesher * cloned = mesher.clone();
      assert_equal(cloned->getUseSimplicesDecomposition(), mesher.getUseSimplicesDecomposition());
      assert(cloned != &mesher);
      delete cloned;
      fullprint << "OK" << std::endl;
    }

    // Test save/load via Study
    fullprint << "--- save/load ---" << std::endl;
    {
      const String studyFile("cdm_test.xml");
      {
        Study study(studyFile);
        study.add("mesher", mesher);
        study.save();
      }
      {
        ConvexDecompositionMesher loadedMesher;
        Study loadedStudy(studyFile);
        loadedStudy.load();
        loadedStudy.fillObject("mesher", loadedMesher);
        assert_equal(loadedMesher.getUseSimplicesDecomposition(), mesher.getUseSimplicesDecomposition());
      }
      std::remove(studyFile.c_str());
      fullprint << "OK" << std::endl;
    }

    // Test invalid mesh dimension != intrinsicDimension (4D mesh with intrinsicDimension=2)
    fullprint << "--- build invalid mesh (dim=4, intrinsicDim=2 via LevelSet) ---" << std::endl;
    {
      const SymbolicFunction f(Description({"x0", "x1", "x2", "x3"}), Description({"(x0^2 + x1^2 + x2^2 + x3^2 + 3)^2 - 16 * (x0^2 + x1^2)"}));
      const LevelSet levelSet(f, LessOrEqual(), 0.0);
      const Mesh mesh = LevelSetMesher(Indices({6, 6, 3, 3})).build(levelSet, Interval(Point({-3.0, -3.0, -1.0, -1.0}), Point({3.0, 3.0, 1.0, 1.0})));
      assert(mesh.isValid());

      // build should work: dim=4, intrinsicDim=2 -> dimension == intrinsicDimension -> simplices decomposition
      Collection<Mesh> decomp = mesher.build(mesh);
      Scalar volumeSum = 0.0;
      for (UnsignedInteger i = 0; i < decomp.getSize(); ++ i)
      {
        assert_equal(ConvexDecompositionMesher::IsConvex(decomp[i]), true);
        volumeSum += decomp[i].getVolume();
      }
      assert_almost_equal(volumeSum, mesh.getVolume());
      fullprint << "Number of components: " << decomp.getSize() << ", volume: " << volumeSum << ", OK" << std::endl;
    }

    fullprint << "All tests passed" << std::endl;
  }
  catch (const TestFailed & ex)
  {
    std::cerr << ex.what() << std::endl;
    return ExitCode::Error;
  }
  catch (const std::exception & ex)
  {
    std::cerr << "Unexpected exception: " << ex.what() << std::endl;
    return ExitCode::Error;
  }

  return ExitCode::Success;
}
