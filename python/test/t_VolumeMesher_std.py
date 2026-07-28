#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing as otm

ot.TESTPREAMBLE()

mesher = otm.VolumeMesher()
print("mesher=", mesher)

assert str(mesher).startswith("class=VolumeMesher")
assert mesher.getApexStrategy() == otm.VolumeMesher.CENTROID

mesher.setApexStrategy(otm.VolumeMesher.FIRST_VERTEX)
assert mesher.getApexStrategy() == otm.VolumeMesher.FIRST_VERTEX
mesher.setApexStrategy(otm.VolumeMesher.CENTROID)
assert mesher.getApexStrategy() == otm.VolumeMesher.CENTROID

hull_mesher = otm.ConvexHullMesher()

# 1. empty mesh -> error
with ott.assert_raises(Exception):
    mesher.build(ot.Mesh())

# 2. tetrahedron surface -> CENTROID
tetra_pts = ot.Sample([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
tetra_surface = hull_mesher.build(tetra_pts)

mesher.setApexStrategy(otm.VolumeMesher.CENTROID)
vol = mesher.build(tetra_surface)
assert vol.getDimension() == 3
assert vol.isValid()
assert vol.isConvex()
ott.assert_almost_equal(vol.getVolume(), 1.0 / 6.0, 1e-3, 1e-3)

# 3. tetrahedron surface -> FIRST_VERTEX
mesher.setApexStrategy(otm.VolumeMesher.FIRST_VERTEX)
vol_first = mesher.build(tetra_surface)
assert vol_first.getDimension() == 3
assert vol_first.isValid()
assert vol.isConvex()
ott.assert_almost_equal(vol_first.getVolume(), 1.0 / 6.0, 1e-3, 1e-3)

# 4. unit cube surface -> CENTROID
cube_corners = ot.IntervalMesher([1] * 3).build(ot.Interval(3)).getVertices()
cube_surface = hull_mesher.build(cube_corners)

mesher.setApexStrategy(otm.VolumeMesher.CENTROID)
cube_vol = mesher.build(cube_surface)
assert cube_vol.getDimension() == 3
assert cube_vol.isValid()
assert cube_vol.isConvex()
ott.assert_almost_equal(cube_vol.getVolume(), 1.0, 1e-3, 1e-3)

# 5. both strategies give same volume
mesher.setApexStrategy(otm.VolumeMesher.CENTROID)
vol_c = mesher.build(tetra_surface)
mesher.setApexStrategy(otm.VolumeMesher.FIRST_VERTEX)
vol_f = mesher.build(tetra_surface)
ott.assert_almost_equal(vol_c.getVolume(), vol_f.getVolume(), 1e-3, 1e-3)

# 6. volume mesh -> surface (via BoundaryMesher) -> volume (round-trip)
volume_mesh = otm.CloudMesher().build(tetra_pts)
surface_from_volume = ot.BoundaryMesher().build(volume_mesh)
surface_from_volume.setIsConvex(volume_mesh.isConvex())  # TODO: drop for OT 1.28
mesher.setApexStrategy(otm.VolumeMesher.CENTROID)
roundtrip_vol = mesher.build(surface_from_volume)
assert roundtrip_vol.getDimension() == 3
assert roundtrip_vol.isValid()
assert roundtrip_vol.isConvex()
ott.assert_almost_equal(roundtrip_vol.getVolume(), 1.0 / 6.0, 1e-3, 1e-3)

# 7. 2D triangle surface -> volume (arbitrary dimension support)
pts2d = ot.Sample([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
surface2d = hull_mesher.build(pts2d)
vol2d = mesher.build(surface2d)
assert vol2d.getDimension() == 2
assert vol2d.isValid()
assert vol2d.isConvex()
ott.assert_almost_equal(vol2d.getVolume(), 0.5, 1e-3, 1e-3)

# 8. reuse mesher
vol1 = mesher.build(tetra_surface)
vol2 = mesher.build(cube_surface)
print("tetra tets=", vol1.getSimplicesNumber())
print("cube  tets=", vol2.getSimplicesNumber())
print("2d tris=", vol2d.getSimplicesNumber())
