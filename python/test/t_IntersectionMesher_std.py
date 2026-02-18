#! /usr/bin/env python

import math
import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.IntersectionMesher()
print("mesher=", mesher)

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

# 3-d example
N = 10
M = 2
dim = 3
xc1 = 0.0
yc1 = 0.0
R1 = 1.0
H1 = 4.0
yc2 = 0.8
zc2 = 0.0
R2 = 0.5
H2 = 4.0
nTheta = 16
cirle1 = []
cirle2 = []
for i in range(nTheta):
    theta = i * 2.0 * math.pi / nTheta
    x1 = R1 * math.cos(theta)
    y1 = R1 * math.sin(theta)
    cirle1.append([x1, y1])
    y2 = R2 * math.cos(theta)
    z2 = R2 * math.sin(theta)
    cirle2.append([y2, z2])
disc1 = otmeshing.PolygonMesher().build(cirle1)
disc2 = otmeshing.PolygonMesher().build(cirle2)

extension1 = ot.Interval([-H1 / 2] * (dim - 2), [H1 / 2] * (dim - 2))
injection1 = [2]  # add z component
cyl1 = otmeshing.Cylinder(disc1, extension1, injection1, M)

extension2 = ot.Interval([-H2 / 2] * (dim - 2), [H2 / 2] * (dim - 2))
injection2 = [0]  # add x component
cyl2 = otmeshing.Cylinder(disc2, extension2, injection2, M)

mesh1 = otmeshing.CloudMesher().build(cyl1.getVertices())
mesh2 = otmeshing.CloudMesher().build(cyl2.getVertices())
# mesh1.exportToVTKFile("mesh1.vtk")
# mesh2.exportToVTKFile("mesh2.vtk")

inter12 = otmeshing.IntersectionMesher().build([mesh1, mesh2])
volume = inter12.getVolume()
print("inter volume=", volume)
ott.assert_almost_equal(volume, 1.462974)
# inter12.exportToVTKFile("inter12.vtk")

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
