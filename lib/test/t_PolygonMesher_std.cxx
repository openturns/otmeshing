#include <iostream>
#include <cmath>

#include "otmeshing/otmeshing.hxx"
#include "otmeshing/PolygonMesher.hxx"

using namespace OT;
using namespace OTMESHING;

int main()
{
  PolygonMesher mesher;
  std::cout << "mesher=" << mesher << std::endl;

  // 2D convex polygon: unit square
  Sample points(4, 2);
  points(0, 0) = 0.0; points(0, 1) = 0.0;
  points(1, 0) = 1.0; points(1, 1) = 0.0;
  points(2, 0) = 1.0; points(2, 1) = 1.0;
  points(3, 0) = 0.0; points(3, 1) = 1.0;

  Mesh mesh(mesher.build(points));
  std::cout << "mesh=" << mesh << std::endl;
  if (mesh.getDimension() != 2)
    throw std::runtime_error("Expected dimension 2");
  if (mesh.getSimplices().getSize() != 2)
    throw std::runtime_error("Expected 2 simplices");
  if (!mesh.isValid())
    throw std::runtime_error("Expected valid mesh");
  if (std::abs(mesh.getVolume() - 1.0) > 1e-14)
    throw std::runtime_error("Expected volume 1.0");

  // 3D coplanar points (z=0)
  Sample points3d(4, 3);
  points3d(0, 0) = 0.0; points3d(0, 1) = 0.0; points3d(0, 2) = 0.0;
  points3d(1, 0) = 1.0; points3d(1, 1) = 0.0; points3d(1, 2) = 0.0;
  points3d(2, 0) = 1.0; points3d(2, 1) = 1.0; points3d(2, 2) = 0.0;
  points3d(3, 0) = 0.0; points3d(3, 1) = 1.0; points3d(3, 2) = 0.0;

  Mesh mesh3d(mesher.build(points3d));
  std::cout << "mesh3d=" << mesh3d << std::endl;
  if (mesh3d.getDimension() != 3)
    throw std::runtime_error("Expected dimension 3");
  if (mesh3d.getSimplices().getSize() != 2)
    throw std::runtime_error("Expected 2 simplices");
  if (!mesh3d.isValid())
    throw std::runtime_error("Expected valid mesh");

  // triangle (minimal case)
  Sample triangle(3, 2);
  triangle(0, 0) = 0.0; triangle(0, 1) = 0.0;
  triangle(1, 0) = 1.0; triangle(1, 1) = 0.0;
  triangle(2, 0) = 0.0; triangle(2, 1) = 1.0;

  Mesh triMesh(mesher.build(triangle));
  std::cout << "triMesh=" << triMesh << std::endl;
  if (triMesh.getSimplices().getSize() != 1)
    throw std::runtime_error("Expected 1 simplex");
  if (!triMesh.isValid())
    throw std::runtime_error("Expected valid mesh");

  std::cout << "ok" << std::endl;
  return 0;
}
