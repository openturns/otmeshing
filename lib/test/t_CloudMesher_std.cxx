#include <iostream>
#include <cmath>

#include <openturns/OT.hxx>
#include <openturns/OTtestcode.hxx>

#include "otmeshing/otmeshing.hxx"

using namespace OT;
using namespace OT::Test;
using namespace OTMESHING;

int main()
{
  // Test default constructor and string representation
  CloudMesher mesher;
  std::cout << "mesher=" << mesher << std::endl;
  assert_equal(mesher.getClassName(), String("CloudMesher"));

  // Test 2D triangulation with BASIC (default)
  {
    Sample points(5, 2);
    points[0] = Point({3.0, 0.0});
    points[1] = Point({2.0, 0.0});
    points[2] = Point({2.0, 0.75});
    points[3] = Point({2.5, 0.75});
    points[4] = Point({3.0, 0.2});
    Mesh mesh = mesher.build(points);
    assert_equal(mesh.getDimension(), 2ul);
    assert_equal(mesh.getSimplices().getSize(), 3ul);
    assert_equal(mesh.isValid(), true);
    assert_equal(mesh.isConvex(), true);
    assert_almost_equal(mesh.getVolume(), 0.6125, 1e-12, 1e-12);
  }

  // Test 2D triangulation with DELAUNAY
  {
    CloudMesher delMesher(CloudMesher::DELAUNAY);
    Sample points(5, 2);
    points[0] = Point({3.0, 0.0});
    points[1] = Point({2.0, 0.0});
    points[2] = Point({2.0, 0.75});
    points[3] = Point({2.5, 0.75});
    points[4] = Point({3.0, 0.2});
    Mesh mesh = delMesher.build(points);
    assert_equal(mesh.getDimension(), 2ul);
    assert_equal(mesh.getSimplices().getSize(), 3ul);
    assert_equal(mesh.isValid(), true);
    assert_equal(mesh.isConvex(), true);
    assert_almost_equal(mesh.getVolume(), 0.6125, 1e-12, 1e-12);
  }

  // Test 1D case: only min and max become vertices
  {
    Sample points(3, 1);
    points[0] = Point({2.5});
    points[1] = Point({1.5});
    points[2] = Point({3.0});
    Mesh mesh = mesher.build(points);
    assert_equal(mesh.getDimension(), 1ul);
    assert_equal(mesh.isValid(), true);
    assert_equal(mesh.isConvex(), true);
    assert_almost_equal(mesh.getVolume(), 1.5, 1e-12, 1e-12);
    assert_equal(mesh.getVertices().getSize(), 2ul);
    assert_equal(mesh.getSimplices().getSize(), 1ul);
    assert_almost_equal(mesh.getVertices()[0][0], 1.5, 1e-12, 1e-12);
    assert_almost_equal(mesh.getVertices()[1][0], 3.0, 1e-12, 1e-12);
  }

  // Test unit hypercube triangulations
  {
    for (UnsignedInteger method = 0; method <= 1; ++method)
    {
      CloudMesher cubeMesher(static_cast<CloudMesher::TriangulationMethod>(method));
      for (UnsignedInteger dim = 1; dim <= 6; ++dim)
      {
        // Generate unit hypercube vertices
        const UnsignedInteger nbVertices = static_cast<UnsignedInteger>(std::pow(2.0, static_cast<double>(dim)));
        Sample vertices(nbVertices, dim);
        for (UnsignedInteger i = 0; i < nbVertices; ++i)
        {
          for (UnsignedInteger j = 0; j < dim; ++j)
          {
            vertices[i][j] = (i >> j) & 1;
          }
        }
        Mesh mesh = cubeMesher.build(vertices);
        assert_equal(mesh.getDimension(), dim);
        assert_equal(mesh.getVertices().getSize(), vertices.getSize());
        assert_equal(mesh.isValid(), true);
        assert_equal(mesh.isConvex(), true);
        assert_almost_equal(mesh.getVolume(), 1.0, 1e-12, 1e-12);
      }
    }
  }

  // Error cases: null dimension
  {
    Sample empty;
    bool caught = false;
    try
    {
      mesher.build(empty);
    }
    catch (const InvalidArgumentException &)
    {
      caught = true;
    }
    assert_equal(caught, true);
  }

  // Error cases: insufficient points for 2D
  {
    Sample tooFew(1, 2);
    tooFew[0] = Point({1.0, 2.0});
    bool caught = false;
    try
    {
      mesher.build(tooFew);
    }
    catch (const InvalidArgumentException &)
    {
      caught = true;
    }
    assert_equal(caught, true);
  }

  // Error cases: 1D with single point
  {
    Sample tooFew(1, 1);
    tooFew[0] = Point({1.0});
    bool caught = false;
    try
    {
      mesher.build(tooFew);
    }
    catch (const InvalidArgumentException &)
    {
      caught = true;
    }
    assert_equal(caught, true);
  }

  std::cout << "All C++ tests passed!" << std::endl;
  return 0;
}
