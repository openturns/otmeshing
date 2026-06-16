//                                               -*- C++ -*-
/**
 *  @brief The test file of class Cylinder for standard methods
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

#include <cstdio>

#include <openturns/OT.hxx>
#include <openturns/OTtestcode.hxx>

#include "otmeshing/Cylinder.hxx"
#include "otmeshing/IntersectionMesher.hxx"

using namespace OT;
using namespace OT::Test;
using namespace OTMESHING;

int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {
    // Test default constructor
    fullprint << "--- Default constructor ---" << std::endl;
    Cylinder cylDefault;
    assert_equal(cylDefault.getDimension(), UnsignedInteger(0));
    fullprint << "Default dimension: " << cylDefault.getDimension() << std::endl;

    // Build a 2D base mesh: square [-1,1]^2 divided into 5x5 cells
    const Indices baseDisc({5, 5});
    const Interval baseInterval(Point({-1.0, -1.0}), Point({1.0, 1.0}));
    const Mesh baseMesh(IntervalMesher(baseDisc).build(baseInterval));
    fullprint << "Base mesh vertices: " << baseMesh.getVerticesNumber() << std::endl;
    fullprint << "Base mesh volume: " << baseMesh.getVolume() << std::endl;

    // Extension: 1D interval [-2, 2]
    const Interval extension(Point({-2.0}), Point({2.0}));
    fullprint << "Extension volume: " << extension.getVolume() << std::endl;

    // Injection: extension maps to dimension 2 (0-indexed)
    const Indices injection({2});

    // Discretization of extension dimension
    const UnsignedInteger M = 3;

    // Test parameter constructor
    fullprint << "--- Parameter constructor ---" << std::endl;
    const Cylinder cylinder(baseMesh, extension, injection, M);
    fullprint << "cylinder: " << cylinder << std::endl;

    // Test getDimension
    fullprint << "Dimension: " << cylinder.getDimension() << std::endl;
    assert_equal(cylinder.getDimension(), UnsignedInteger(3));

    // Test getBase
    fullprint << "--- getBase ---" << std::endl;
    const Mesh retrievedBase = cylinder.getBase();
    assert_almost_equal(retrievedBase.getVertices(), baseMesh.getVertices());
    assert_equal(retrievedBase.getSimplices(), baseMesh.getSimplices());

    // Test getExtension
    fullprint << "--- getExtension ---" << std::endl;
    const Interval retrievedExt = cylinder.getExtension();
    assert_equal(retrievedExt.getLowerBound()[0], extension.getLowerBound()[0]);
    assert_equal(retrievedExt.getUpperBound()[0], extension.getUpperBound()[0]);

    // Test getInjection
    fullprint << "--- getInjection ---" << std::endl;
    const Indices retrievedInj = cylinder.getInjection();
    assert_equal(retrievedInj.getSize(), injection.getSize());
    assert_equal(retrievedInj[0], injection[0]);

    // Test getDiscretization
    fullprint << "--- getDiscretization ---" << std::endl;
    assert_equal(cylinder.getDiscretization(), M);

    // Test getVertices
    fullprint << "--- getVertices ---" << std::endl;
    const Sample vertices = cylinder.getVertices();
    const UnsignedInteger extVerticesNumber = IntervalMesher(Indices({M})).build(extension).getVertices().getSize();
    fullprint << "Base vertices: " << baseMesh.getVerticesNumber()
              << ", ext vertices: " << extVerticesNumber << std::endl;
    assert_equal(vertices.getSize(), baseMesh.getVerticesNumber() * extVerticesNumber);
    assert_equal(vertices.getDimension(), UnsignedInteger(3));

    // Test getBoundingBox
    fullprint << "--- getBoundingBox ---" << std::endl;
    const Interval bbox = cylinder.getBoundingBox();
    fullprint << "BBox: " << bbox << std::endl;
    assert_almost_equal(bbox.getLowerBound()[0], -1.0);
    assert_almost_equal(bbox.getLowerBound()[1], -1.0);
    assert_almost_equal(bbox.getLowerBound()[2], -2.0);
    assert_almost_equal(bbox.getUpperBound()[0],  1.0);
    assert_almost_equal(bbox.getUpperBound()[1],  1.0);
    assert_almost_equal(bbox.getUpperBound()[2],  2.0);

    // Test getVolume: base area (4) * extension length (4) = 16
    fullprint << "--- getVolume ---" << std::endl;
    const Scalar expectedVolume = baseMesh.getVolume() * extension.getVolume();
    fullprint << "Expected volume: " << expectedVolume
              << ", cylinder volume: " << cylinder.getVolume() << std::endl;
    assert_almost_equal(cylinder.getVolume(), expectedVolume);

    // Test isConvex (a square-based cylinder is convex)
    fullprint << "--- isConvex ---" << std::endl;
    assert(cylinder.isConvex());
    fullprint << "Convex: true" << std::endl;

    // Test computeMesh
    fullprint << "--- computeMesh ---" << std::endl;
    const Mesh computedMesh = cylinder.computeMesh();
    fullprint << "Computed mesh vertices: " << computedMesh.getVerticesNumber() << std::endl;
    assert(computedMesh.getVerticesNumber() > 0);

    // Test clone
    fullprint << "--- clone ---" << std::endl;
    Cylinder * clonedCylinder = cylinder.clone();
    assert_equal(clonedCylinder->getDimension(), cylinder.getDimension());
    assert_almost_equal(clonedCylinder->getVolume(), cylinder.getVolume());
    assert_almost_equal(clonedCylinder->getBoundingBox().getLowerBound(),
                        cylinder.getBoundingBox().getLowerBound());
    fullprint << "Clone OK" << std::endl;
    delete clonedCylinder;

    // Test string representation
    fullprint << "--- string representation ---" << std::endl;
    {
      const String repr = cylinder.__repr__();
      fullprint << "__repr__: " << repr << std::endl;
      assert_equal(repr, String("class=Cylinder"));
    }
    {
      const String str = cylinder.__str__();
      fullprint << "__str__: " << str << std::endl;
      assert_equal(str, String("class=Cylinder"));
    }

    // Test save/load via Study
    fullprint << "--- save/load ---" << std::endl;
    const String studyFile("cylinder_test.xml");
    {
      Study study(studyFile);
      study.add("cylinder", cylinder);
      study.save();
    }
    {
      Cylinder loadedCylinder;
      Study loadedStudy(studyFile);
      loadedStudy.load();
      loadedStudy.fillObject("cylinder", loadedCylinder);
      assert_equal(loadedCylinder.getDimension(), cylinder.getDimension());
      assert_almost_equal(loadedCylinder.getVolume(), cylinder.getVolume());
      assert_almost_equal(loadedCylinder.getBoundingBox().getLowerBound(),
                          cylinder.getBoundingBox().getLowerBound());
      assert_almost_equal(loadedCylinder.getBoundingBox().getUpperBound(),
                          cylinder.getBoundingBox().getUpperBound());
      assert_equal(loadedCylinder.getDiscretization(), cylinder.getDiscretization());
    }
    std::remove(studyFile.c_str());
    fullprint << "Save/load OK" << std::endl;

    // Test error cases: invalid injection size (size != extension dimension)
    fullprint << "--- error cases ---" << std::endl;
    {
      const Indices badInjection({2, 3});
      bool caught = false;
      try
      {
        const Cylinder badCyl(baseMesh, extension, badInjection, M);
      }
      catch (const InvalidArgumentException &)
      {
        caught = true;
      }
      assert(caught);
      (void)caught;
      fullprint << "Caught invalid injection size" << std::endl;
    }

    {
      const Indices badInjection2({5});
      bool caught = false;
      try
      {
        const Cylinder badCyl2(baseMesh, extension, badInjection2, M);
      }
      catch (const InvalidArgumentException &)
      {
        caught = true;
      }
      assert(caught);
      (void)caught;
      fullprint << "Caught invalid injection index" << std::endl;
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
