#! /usr/bin/env python

import math
import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.IntersectionMesher()
print("mesher=", mesher)
print(mesher.getRecompress())

# intersection of two cubes
for compression in [False, True]:
    mesher.setRecompress(compression)
    for dim in range(2, 6):
        mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [3.0] * dim))
        mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([1.0] * dim, [4.0] * dim))
        intersection = mesher.build([mesh1, mesh2])
        volume = intersection.getVolume()
        print(f"{dim=} {compression=} intersection={intersection} {volume=:.3g}")
        ott.assert_almost_equal(volume, 2.0**dim)

        if (dim == 3) and compression:
            bmesh = ot.BoundaryMesher().build(intersection)
            print(bmesh)
            # bmesh.exportToVTKFile("/tmp/boundary.vtk")

# empty/self intersection
dim = 3
mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval([0.0] * dim, [1.0] * dim))
mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([3.0] * dim, [4.0] * dim))
intersection = mesher.build([mesh1, mesh2])
volume = intersection.getVolume()
print(f"{dim=} {compression=} intersection={intersection} {volume=:.3g}")
ott.assert_almost_equal(volume, 0.0)
intersection = mesher.build([])
volume = intersection.getVolume()
print(f"{dim=} {compression=} intersection={intersection} {volume=:.3g}")
ott.assert_almost_equal(volume, 0.0)
intersection = mesher.build([mesh1])
volume = intersection.getVolume()
print(f"{dim=} {compression=} intersection={intersection} {volume=:.3g}")
assert volume == mesh1.getVolume()
intersection = mesher.build([mesh1, mesh1])
volume = intersection.getVolume()
print(f"{dim=} {compression=} intersection={intersection} {volume=:.3g}")
ott.assert_almost_equal(volume, mesh1.getVolume())

# 3-d cylinder intersection: Steinmetz solid
M = 2
dim = 3
R = 2.0
H = 2.5 * R
nTheta = 32
cirle1 = []
cirle2 = []
star3 = []
for i in range(nTheta):
    theta = i * 2.0 * math.pi / nTheta
    x1 = R * math.cos(theta)
    y1 = R * math.sin(theta)
    cirle1.append([x1, y1])
    y2 = R * math.cos(theta)
    z2 = R * math.sin(theta)
    cirle2.append([y2, z2])
    star3.append([x1, y1] if i % 2 == 0 else [x1 * 1.5, y1 * 1.5])
disc1 = otmeshing.PolygonMesher().build(cirle1)
disc2 = otmeshing.PolygonMesher().build(cirle2)
disc3 = otmeshing.PolygonMesher().build(star3)

extension1 = ot.Interval([-H / 2] * (dim - 2), [H / 2] * (dim - 2))
injection1 = [2]  # add z component
cyl1 = otmeshing.Cylinder(disc1, extension1, injection1, M)
# cylinder3 is non-convex (star-shaped base)
cyl3 = otmeshing.Cylinder(disc3, extension1, injection1, M)

extension2 = ot.Interval([-H / 2] * (dim - 2), [H / 2] * (dim - 2))
injection2 = [0]  # add x component
cyl2 = otmeshing.Cylinder(disc2, extension2, injection2, M)

mesh1 = otmeshing.CloudMesher().build(cyl1.getVertices())
mesh2 = otmeshing.CloudMesher().build(cyl2.getVertices())
# mesh1.exportToVTKFile("mesh1.vtk")
# mesh2.exportToVTKFile("mesh2.vtk")
volume_ref = 16.0 / 3.0 * R**3

# compute intersection with buildConvexSample
assert cyl1.isConvex() and cyl2.isConvex()
interSample = otmeshing.IntersectionMesher().buildConvexSample([cyl1.getVertices(), cyl2.getVertices()])
inter12 = otmeshing.CloudMesher().build(interSample)
volume = inter12.getVolume()
print("inter(sample) volume=", volume)
ott.assert_almost_equal(volume, volume_ref, 1e-2)

# compute intersection with buildCylinder
inter12 = otmeshing.IntersectionMesher().buildCylinder([cyl1, cyl2])
volume = inter12.getVolume()
print("inter(cylinder) volume=", volume)
ott.assert_almost_equal(volume, volume_ref, 1e-2)
# inter12.exportToVTKFile("inter12.vtk")

# compute non-convex intersection with buildCylinder
inter32 = otmeshing.IntersectionMesher().buildCylinder([cyl3, cyl2])
volume = inter32.getVolume()
print("inter(cylinder, non-convex) volume=", volume)
ott.assert_almost_equal(volume, 53.4976)

# convex intersection
for dim in range(2, 6):
    for N in range(1, 25, 5):
        discr = [N, N] + [1] * (dim - 2)
        mesh1 = ot.IntervalMesher(discr).build(ot.Interval([0.0] * dim, [3.0] * dim))
        mesh2 = ot.IntervalMesher(discr).build(ot.Interval([1.0] * dim, [4.0] * dim))
        intersection = mesher.buildConvex([mesh1, mesh2])
        volume = intersection.getVolume()
        print(f"{dim=} {N=} intersection={intersection} {volume=:.3g}")
        ott.assert_almost_equal(volume, 2.0**dim)

# empty/self convex intersection
dim = 3
discr = [2] * dim
mesh1 = ot.IntervalMesher(discr).build(ot.Interval([0.0] * dim, [1.0] * dim))
mesh2 = ot.IntervalMesher(discr).build(ot.Interval([3.0] * dim, [4.0] * dim))
intersection = mesher.buildConvex([mesh1, mesh2])
volume = intersection.getVolume()
print(f"{dim=} intersection={intersection} {volume=:.3g}")
ott.assert_almost_equal(volume, 0.0)
intersection = mesher.buildConvex([])
volume = intersection.getVolume()
print(f"{dim=} intersection={intersection} {volume=:.3g}")
ott.assert_almost_equal(volume, 0.0)
intersection = mesher.buildConvex([mesh1])
volume = intersection.getVolume()
print(f"{dim=} intersection={intersection} {volume=:.3g}")
assert volume == mesh1.getVolume()
intersection = mesher.buildConvex([mesh1, mesh1])
volume = intersection.getVolume()
print(f"{dim=} intersection={intersection} {volume=:.3g}")
ott.assert_almost_equal(volume, mesh1.getVolume())
