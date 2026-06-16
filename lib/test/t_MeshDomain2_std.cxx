#include <iostream>
#include <cmath>

#include "otmeshing/otmeshing.hxx"
#include "otmeshing/MeshDomain2.hxx"

using namespace OT;
using namespace OTMESHING;

int main()
{
  // 2D unit square
  Sample vertices(4, 2);
  vertices(0, 0) = 0.0; vertices(0, 1) = 0.0;
  vertices(1, 0) = 1.0; vertices(1, 1) = 0.0;
  vertices(2, 0) = 1.0; vertices(2, 1) = 1.0;
  vertices(3, 0) = 0.0; vertices(3, 1) = 1.0;
  IndicesCollection simplices(2, 3);
  simplices(0, 0) = 0; simplices(0, 1) = 1; simplices(0, 2) = 2;
  simplices(1, 0) = 0; simplices(1, 1) = 2; simplices(1, 2) = 3;
  Mesh mesh2d(vertices, simplices);

  MeshDomain2 domain2d(mesh2d);
  std::cout << "domain2d=" << domain2d << std::endl;

  // inside point
  Point pInside2d({0.2, 0.2});
  Scalar d = domain2d.computeDistance(pInside2d);
  std::cout << "inside 2d distance=" << d << std::endl;
  if (std::abs(d + 0.2) > 1e-14)
    throw std::runtime_error("Expected distance -0.2");

  // outside point
  Point pOutside2d({1.2, 1.2});
  d = domain2d.computeDistance(pOutside2d);
  std::cout << "outside 2d distance=" << d << std::endl;
  if (std::abs(d - 0.2 * std::sqrt(2.0)) > 1e-14)
    throw std::runtime_error("Expected distance 0.2*sqrt(2)");

  // 3D unit cube
  Sample vertices3d(8, 3);
  for (UnsignedInteger i = 0; i < 8; ++ i)
    for (UnsignedInteger j = 0; j < 3; ++ j)
      vertices3d(i, j) = (i >> j) & 1;
  IndicesCollection simplices3d(6, 4);
  // tetrahedralization of a cube as 6 tets
  simplices3d(0, 0) = 0; simplices3d(0, 1) = 1; simplices3d(0, 2) = 3; simplices3d(0, 3) = 7;
  simplices3d(1, 0) = 0; simplices3d(1, 1) = 1; simplices3d(1, 2) = 5; simplices3d(1, 3) = 7;
  simplices3d(2, 0) = 0; simplices3d(2, 1) = 2; simplices3d(2, 2) = 3; simplices3d(2, 3) = 7;
  simplices3d(3, 0) = 0; simplices3d(3, 1) = 2; simplices3d(3, 2) = 6; simplices3d(3, 3) = 7;
  simplices3d(4, 0) = 0; simplices3d(4, 1) = 4; simplices3d(4, 2) = 5; simplices3d(4, 3) = 7;
  simplices3d(5, 0) = 0; simplices3d(5, 1) = 4; simplices3d(5, 2) = 6; simplices3d(5, 3) = 7;
  Mesh mesh3d(vertices3d, simplices3d);

  MeshDomain2 domain3d(mesh3d);
  std::cout << "domain3d=" << domain3d << std::endl;

  // inside point
  Point pInside3d({0.2, 0.2, 0.2});
  d = domain3d.computeDistance(pInside3d);
  std::cout << "inside 3d distance=" << d << std::endl;
  if (std::abs(d + 0.2) > 1e-14)
    throw std::runtime_error("Expected distance -0.2");

  // outside point
  Point pOutside3d({1.2, 1.2, 1.2});
  d = domain3d.computeDistance(pOutside3d);
  std::cout << "outside 3d distance=" << d << std::endl;
  if (std::abs(d - 0.2 * std::sqrt(3.0)) > 1e-14)
    throw std::runtime_error("Expected distance 0.2*sqrt(3)");

  std::cout << "ok" << std::endl;
  return 0;
}
