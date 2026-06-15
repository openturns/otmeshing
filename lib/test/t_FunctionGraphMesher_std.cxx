//                                               -*- C++ -*-
/**
 *  @brief Test for FunctionGraphMesher
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

#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"

#include "otmeshing/FunctionGraphMesher.hxx"

using namespace OT;
using namespace OT::Test;
using namespace OTMESHING;

int main()
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {

    // 1. Default constructor + string representation
    fullprint << "Default constructor" << std::endl;
    FunctionGraphMesher mesher;
    fullprint << "mesher = " << mesher << std::endl;

    // 2. 1D linear function -- exact volume
    fullprint << "1D linear function" << std::endl;
    SymbolicFunction f(Description({"x"}), Description({"x"}));
    f.setName("linear");
    f.setInputDescription({"x"});
    f.setOutputDescription({"y"});

    FunctionGraphMesher mesher1(Interval({0.0}, {1.0}), Indices({5}));
    fullprint << "mesher = " << mesher1 << std::endl;

    // subGraph=True, outputIndex=1 (output last)
    Mesh mesh(mesher1.build(f, 1, -1.0, 2.0, 1, true));
    fullprint << "mesh (subGraph=T, outIdx=1): " << mesh << std::endl;
    assert_equal(mesh.getDimension(), UnsignedInteger(2));
    assert_equal(mesh.getVerticesNumber(), UnsignedInteger(12));  // (5+1)*(1+1)
    assert_equal(mesh.getSimplicesNumber(), UnsignedInteger(10)); // 5*1*2
    assert_almost_equal(mesh.getVolume(), 1.5);  // int(x+1)dx = 3/2
    assert_equal(mesh.getDescription(), Description({"x", "y"}));

    // subGraph=False
    Mesh mesh2(mesher1.build(f, 1, -1.0, 2.0, 1, false));
    fullprint << "mesh (subGraph=F, outIdx=1): vol=" << mesh2.getVolume() << std::endl;
    assert_almost_equal(mesh2.getVolume(), 1.5);  // int(2-x)dx = 3/2

    // outputIndex=0 (output first)
    Mesh mesh3(mesher1.build(f, 0, -1.0, 2.0, 1, true));
    fullprint << "mesh (subGraph=T, outIdx=0): dim=" << mesh3.getDimension() << " vol=" << mesh3.getVolume() << std::endl;
    assert_equal(mesh3.getDimension(), UnsignedInteger(2));
    assert_almost_equal(mesh3.getVolume(), 1.5);
    assert_equal(mesh3.getDescription(), Description({"y", "x"}));

    // outputDiscretization > 1
    Mesh mesh4(mesher1.build(f, 1, -1.0, 2.0, 3, true));
    fullprint << "mesh (subGraph=T, outDisc=3): vert=" << mesh4.getVerticesNumber() << " simp=" << mesh4.getSimplicesNumber() << " vol=" << mesh4.getVolume() << std::endl;
    assert_equal(mesh4.getVerticesNumber(), UnsignedInteger(24));   // (5+1)*(3+1)
    assert_equal(mesh4.getSimplicesNumber(), UnsignedInteger(30));  // 5*3*2
    assert_almost_equal(mesh4.getVolume(), 1.5);

    // 3. 1D constant function
    fullprint << "1D constant function" << std::endl;
    SymbolicFunction f_zero(Description({"x"}), Description({"0.0"}));
    FunctionGraphMesher mesher_cst(Interval({0.0}, {1.0}), Indices({1}));
    mesh = mesher_cst.build(f_zero, 1, -1.0, 1.0, 1, true);
    assert_almost_equal(mesh.getVolume(), 1.0);  // int(0-(-1))dx = 1
    mesh = mesher_cst.build(f_zero, 1, -1.0, 1.0, 1, false);
    assert_almost_equal(mesh.getVolume(), 1.0);  // int(1-0)dx = 1

    // 4. 2D linear function -- exact volume
    fullprint << "2D linear function" << std::endl;
    SymbolicFunction f2d(Description({"x0", "x1"}), Description({"x0 + x1"}));
    f2d.setInputDescription({"x", "y"});
    f2d.setOutputDescription({"z"});
    FunctionGraphMesher mesher3(Interval({0.0, 0.0}, {1.0, 1.0}), Indices({1, 1}));
    mesh = mesher3.build(f2d, 2, -1.0, 3.0, 1, true);
    fullprint << "mesh (subGraph=T, outIdx=2): dim=" << mesh.getDimension() << " vol=" << mesh.getVolume() << std::endl;
    assert_equal(mesh.getDimension(), UnsignedInteger(3));
    assert_almost_equal(mesh.getVolume(), 2.0);  // int(x0+x1+1) = 2

    // subGraph=False
    mesh = mesher3.build(f2d, 2, -1.0, 3.0, 1, false);
    fullprint << "  supergraph volume: " << mesh.getVolume() << std::endl;
    assert_almost_equal(mesh.getVolume(), 2.0);  // int(3-x0-x1) = 2

    // outputIndex at other positions
    for (UnsignedInteger idx = 0; idx < 2; ++idx)
    {
      mesh = mesher3.build(f2d, idx, -1.0, 3.0, 1, true);
      fullprint << "  outputIndex=" << idx << " volume: " << mesh.getVolume() << std::endl;
      assert_almost_equal(mesh.getVolume(), 2.0);
    }

    // outputDiscretization > 1
    mesh = mesher3.build(f2d, 2, -1.0, 3.0, 3, true);
    fullprint << "  outDisc=3: vert=" << mesh.getVerticesNumber() << " simp=" << mesh.getSimplicesNumber() << " vol=" << mesh.getVolume() << std::endl;
    assert_equal(mesh.getVerticesNumber(), UnsignedInteger(16));   // 2*2*4
    assert_equal(mesh.getSimplicesNumber(), UnsignedInteger(18));  // 1*1*3*6
    assert_almost_equal(mesh.getVolume(), 2.0);

    // 5. 2D non-linear function (original test)
    fullprint << "2D non-linear function" << std::endl;
    Point a({-4.0, -4.0, -4.0});
    Point b({4.0, 4.0, 4.0});
    SymbolicFunction f_orig(Description({"x0", "x1"}), Description({"cos(pi_*x0)*sin(pi_*x1)^2"}));
    f_orig.setName("Paraboloid");
    f_orig.setInputDescription(Description({"$x_0$", "$x_1$"}));
    f_orig.setOutputDescription(Description({"$x_2$"}));

    FunctionGraphMesher mesherOrig(Interval({a[0], a[1]}, {b[0], b[1]}), Indices({100, 100}));
    const UnsignedInteger outputIndex = 2;
    const UnsignedInteger outputDiscretization = 1;
    const Bool subGraph = true;

    mesh = mesherOrig.build(f_orig, outputIndex, a[2], b[2], outputDiscretization, subGraph);
    fullprint << mesh << std::endl;
    assert_equal(mesh.getDimension(), UnsignedInteger(3));
    assert_equal(mesh.getName(), String("Paraboloid"));
    assert_equal(mesh.getVerticesNumber(), UnsignedInteger(20402));
    assert_equal(mesh.getSimplicesNumber(), UnsignedInteger(60000));
    assert_almost_equal(mesh.getVolume(), 256.0);
    assert_equal(mesh.getDescription(), Description({"$x_0$", "$x_1$", "$x_2$"}));

    // subGraph=False -- supergraph volume equals subgraph volume (symmetric f)
    mesh = mesherOrig.build(f_orig, outputIndex, a[2], b[2], outputDiscretization, false);
    fullprint << "  supergraph volume: " << mesh.getVolume() << std::endl;
    assert_almost_equal(mesh.getVolume(), 256.0);

    // outputDiscretization > 1
    mesh = mesherOrig.build(f_orig, outputIndex, a[2], b[2], 2, true);
    fullprint << "  outDisc=2: vert=" << mesh.getVerticesNumber() << " simp=" << mesh.getSimplicesNumber() << " vol=" << mesh.getVolume() << std::endl;
    assert_equal(mesh.getVerticesNumber(), UnsignedInteger(30603));
    assert_equal(mesh.getSimplicesNumber(), UnsignedInteger(120000));
    assert_almost_equal(mesh.getVolume(), 256.0);

    // outputIndex at 0
    mesh = mesherOrig.build(f_orig, 0, a[2], b[2], outputDiscretization, true);
    fullprint << "  outputIndex=0: vol=" << mesh.getVolume() << std::endl;
    assert_almost_equal(mesh.getVolume(), 256.0);

    // 6. Exception tests
    fullprint << "Exception tests" << std::endl;

    // Wrong input dimension
    try
    {
      mesherOrig.build(SymbolicFunction(Description({"x0", "x1", "x2"}), Description({"x0"})), outputIndex, a[2], b[2], outputDiscretization, subGraph);
      throw TestFailed("Should have thrown for wrong input dim");
    }
    catch (const InvalidArgumentException &)
    {
      fullprint << "OK: wrong input dim" << std::endl;
    }

    // Wrong output dimension
    try
    {
      mesherOrig.build(SymbolicFunction(Description({"x0", "x1"}), Description({"x0", "x1"})), outputIndex, a[2], b[2], outputDiscretization, subGraph);
      throw TestFailed("Should have thrown for wrong output dim");
    }
    catch (const InvalidArgumentException &)
    {
      fullprint << "OK: wrong output dim" << std::endl;
    }

    // minOutput >= maxOutput
    try
    {
      mesherOrig.build(f_orig, outputIndex, 0.0, 0.0, outputDiscretization, subGraph);
      throw TestFailed("Should have thrown for equal min/max");
    }
    catch (const InvalidArgumentException &)
    {
      fullprint << "OK: equal min/max" << std::endl;
    }

    // outputIndex > inputDimension
    try
    {
      mesherOrig.build(f_orig, 5, a[2], b[2], outputDiscretization, subGraph);
      throw TestFailed("Should have thrown for large outputIndex");
    }
    catch (const InvalidArgumentException &)
    {
      fullprint << "OK: large outputIndex" << std::endl;
    }

    // outputDiscretization == 0
    try
    {
      mesherOrig.build(f_orig, outputIndex, a[2], b[2], 0, subGraph);
      throw TestFailed("Should have thrown for zero outputDisc");
    }
    catch (const InvalidArgumentException &)
    {
      fullprint << "OK: zero outputDisc" << std::endl;
    }

  }
  catch (TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }

  return ExitCode::Success;
}
