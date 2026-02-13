#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing as otm
import time

ot.TESTPREAMBLE()

N = 20
M = 2
dim = 3
print("Build cylinder 1")
xc1 = 0.0
yc1 = 0.0
R1 = 1.0
H1 = 4.0
f1 = ot.SymbolicFunction(
    ["x", "y"], [f"(x-({xc1}))^2+(y-({yc1}))^2"]
)
levelSet1 = ot.LevelSet(f1, ot.LessOrEqual(), R1**2)
base1 = ot.LevelSetMesher([N] * 2).build(
    levelSet1, ot.Interval([xc1 - R1, yc1 - R1], [xc1 + R1, yc1 + R1])
)
extension1 = ot.Interval([-H1 / 2] * (dim - 2), [H1 / 2] * (dim - 2))
injection1 = list(range(2, dim))
cyl1 = otm.Cylinder(base1, extension1, injection1, M)

print("Build cylinder 2")
yc2 = 0.8
zc2 = 0.0
R2 = 0.5
H2 = 4.0
f2 = ot.SymbolicFunction(
    ["y", "z"], [f"(y-({yc2}))^2+(z-({zc2}))^2"]
)
levelSet2 = ot.LevelSet(f2, ot.LessOrEqual(), R2**2)
base2 = ot.LevelSetMesher([N] * 2).build(
    levelSet2, ot.Interval([yc2 - R2, zc2 - R2], [yc2 + R2, zc2 + R2])
)
extension2 = ot.Interval([-H2 / 2] * (dim - 2), [H2 / 2] * (dim - 2))
injection2 = [0] + list(range(3, dim))
cyl2 = otm.Cylinder(base2, extension2, injection2, M)

print("Mesh cylinder 1")
method = otm.CloudMesher.BASIC
mesh1 = otm.CloudMesher(method).build(cyl1.getVertices())
if dim <= 3:
    mesh1.exportToVTKFile("mesh1.vtk")
print("Mesh cylinder 2")
mesh2 = otm.CloudMesher(method).build(cyl2.getVertices())
if dim <= 3:
    mesh2.exportToVTKFile("mesh2.vtk")

print("Intersect cylinders")
t0 = time.time()
inter12 = otm.IntersectionMesher().buildConvex([mesh1, mesh2])
t1 = time.time()
print("Convex, t=", t1 - t0, "s", "volume=", inter12.getVolume())
ott.assert_almost_equal(inter12.getVolume(), 0.778142671)
if dim <= 3:
    inter12.exportToVTKFile("inter12_convex.vtk")

quit()
t0 = time.time()
algo = otm.IntersectionMesher()
inter12 = algo.build(mesh1, mesh2)
t1 = time.time()
print("Convex, t=", t1 - t0, "s", "volume=", inter12.getVolume())
if dim <= 3:
    inter12.exportToVTKFile("inter12_general.vtk")
