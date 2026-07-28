//                                               -*- C++ -*-
/**
 *  @brief The test file of class VolumeMesher
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
#include "otmeshing/VolumeMesher.hxx"
#include "otmeshing/ConvexHullMesher.hxx"
#include "otmeshing/CloudMesher.hxx"
#include "openturns/OTtestcode.hxx"
#include "openturns/BoundaryMesher.hxx"
#include "openturns/IntervalMesher.hxx"
#include "openturns/Interval.hxx"

using namespace OT;
using namespace OTMESHING;
using namespace OT::Test;

int main()
{
  TESTPREAMBLE;

  VolumeMesher mesher;
  std::cout << mesher << std::endl;
  assert_equal(mesher.getApexStrategy(), VolumeMesher::CENTROID);

  ConvexHullMesher hull;

  // 1. empty mesh -> error
  try
  {
    mesher.build(Mesh());
    throw TestFailed(OSS() << "Expected InvalidArgumentException");
  }
  catch (const InvalidArgumentException &)
  {
    std::cout << "Empty mesh correctly raised" << std::endl;
  }

  // 2. tetrahedron surface -> CENTROID
  {
    Sample tetra_pts(0, 3);
    tetra_pts.add(Point({0.0, 0.0, 0.0}));
    tetra_pts.add(Point({1.0, 0.0, 0.0}));
    tetra_pts.add(Point({0.0, 1.0, 0.0}));
    tetra_pts.add(Point({0.0, 0.0, 1.0}));
    const Mesh tetra_surface(hull.build(tetra_pts));

    mesher.setApexStrategy(VolumeMesher::CENTROID);
    const Mesh vol(mesher.build(tetra_surface));
    std::cout << "3D tetra centroid: " << vol << std::endl;
    assert_equal(vol.getDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(vol.getIntrinsicDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(vol.getVerticesNumber(), static_cast<UnsignedInteger>(5));
    assert_equal(vol.getSimplicesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(vol.getVolume(), 1.0 / 6.0);
    assert_equal(vol.isValid(), true);
  }

  // 3. tetrahedron surface -> FIRST_VERTEX
  {
    Sample tetra_pts(0, 3);
    tetra_pts.add(Point({0.0, 0.0, 0.0}));
    tetra_pts.add(Point({1.0, 0.0, 0.0}));
    tetra_pts.add(Point({0.0, 1.0, 0.0}));
    tetra_pts.add(Point({0.0, 0.0, 1.0}));
    const Mesh tetra_surface(hull.build(tetra_pts));

    mesher.setApexStrategy(VolumeMesher::FIRST_VERTEX);
    const Mesh vol(mesher.build(tetra_surface));
    std::cout << "3D tetra first vertex: " << vol << std::endl;
    assert_equal(vol.getDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(vol.getIntrinsicDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(vol.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_equal(vol.getSimplicesNumber(), static_cast<UnsignedInteger>(1));
    assert_almost_equal(vol.getVolume(), 1.0 / 6.0);
    assert_equal(vol.isValid(), true);
  }

  // 4. cube surface -> CENTROID
  {
    const Interval cube_interval(3);
    const Sample cube_corners(
      IntervalMesher(Indices(3, 1)).build(cube_interval).getVertices());
    const Mesh cube_surface(hull.build(cube_corners));

    mesher.setApexStrategy(VolumeMesher::CENTROID);
    const Mesh cube_vol(mesher.build(cube_surface));
    std::cout << "cube: " << cube_vol << std::endl;
    assert_equal(cube_vol.getDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(cube_vol.getIntrinsicDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(cube_vol.getSimplicesNumber(), static_cast<UnsignedInteger>(12));
    assert_equal(cube_vol.getVerticesNumber(), static_cast<UnsignedInteger>(9));
    assert_almost_equal(cube_vol.getVolume(), 1.0);
    assert_equal(cube_vol.isValid(), true);
  }

  // 5. both strategies give same volume
  {
    Sample tetra_pts(0, 3);
    tetra_pts.add(Point({0.0, 0.0, 0.0}));
    tetra_pts.add(Point({1.0, 0.0, 0.0}));
    tetra_pts.add(Point({0.0, 1.0, 0.0}));
    tetra_pts.add(Point({0.0, 0.0, 1.0}));
    const Mesh tetra_surface(hull.build(tetra_pts));

    mesher.setApexStrategy(VolumeMesher::CENTROID);
    const Mesh vol_c(mesher.build(tetra_surface));
    mesher.setApexStrategy(VolumeMesher::FIRST_VERTEX);
    const Mesh vol_f(mesher.build(tetra_surface));
    assert_almost_equal(vol_c.getVolume(), vol_f.getVolume());
  }

  // 6. volume mesh -> surface (BoundaryMesher) -> volume (round-trip)
  {
    Sample tetra_pts(0, 3);
    tetra_pts.add(Point({0.0, 0.0, 0.0}));
    tetra_pts.add(Point({1.0, 0.0, 0.0}));
    tetra_pts.add(Point({0.0, 1.0, 0.0}));
    tetra_pts.add(Point({0.0, 0.0, 1.0}));
    const Mesh volume_mesh(CloudMesher().build(tetra_pts));
    Mesh surface_from_volume(BoundaryMesher().build(volume_mesh));
    surface_from_volume.setIsConvex(volume_mesh.isConvex());  // TODO: drop for OT 1.28
    mesher.setApexStrategy(VolumeMesher::CENTROID);
    const Mesh roundtrip_vol(mesher.build(surface_from_volume));
    std::cout << "round-trip: " << roundtrip_vol << std::endl;
    assert_equal(roundtrip_vol.getDimension(), static_cast<UnsignedInteger>(3));
    assert_equal(roundtrip_vol.getIntrinsicDimension(), static_cast<UnsignedInteger>(3));
    assert_almost_equal(roundtrip_vol.getVolume(), 1.0 / 6.0);
    assert_equal(roundtrip_vol.isValid(), true);
  }

  // 7. 2D triangle surface -> volume (arbitrary dimension)
  {
    Sample pts2d(0, 2);
    pts2d.add(Point({0.0, 0.0}));
    pts2d.add(Point({1.0, 0.0}));
    pts2d.add(Point({0.0, 1.0}));
    const Mesh surface2d(hull.build(pts2d));

    mesher.setApexStrategy(VolumeMesher::CENTROID);
    const Mesh vol2d(mesher.build(surface2d));
    std::cout << "2D: " << vol2d << std::endl;
    assert_equal(vol2d.getDimension(), static_cast<UnsignedInteger>(2));
    assert_equal(vol2d.getIntrinsicDimension(), static_cast<UnsignedInteger>(2));
    assert_equal(vol2d.getSimplicesNumber(), static_cast<UnsignedInteger>(3));
    assert_equal(vol2d.getVerticesNumber(), static_cast<UnsignedInteger>(4));
    assert_almost_equal(vol2d.getVolume(), 0.5);
    assert_equal(vol2d.isValid(), true);
  }

  // 8. reuse mesher
  {
    Sample tetra_pts(0, 3);
    tetra_pts.add(Point({0.0, 0.0, 0.0}));
    tetra_pts.add(Point({1.0, 0.0, 0.0}));
    tetra_pts.add(Point({0.0, 1.0, 0.0}));
    tetra_pts.add(Point({0.0, 0.0, 1.0}));
    const Mesh surface(hull.build(tetra_pts));

    mesher.setApexStrategy(VolumeMesher::CENTROID);
    const Mesh vol1(mesher.build(surface));
    const Mesh vol2(mesher.build(surface));
    assert_equal(vol1.getVolume(), vol2.getVolume());
  }

  return 0;
}
