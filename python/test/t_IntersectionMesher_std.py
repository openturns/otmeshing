#! /usr/bin/env python

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

# 3-d example
N = 10
f1 = ot.SymbolicFunction(["x", "y", "z"], ["x^4-y^4+2*(x^2+z^2)-16"])
level1 = ot.LevelSet(f1, ot.LessOrEqual(), 0.0)
mesh1 = ot.LevelSetMesher([N] * 3).build(level1, ot.Interval([-5.0] * 3, [5.0] * 3))
# mesh1.exportToVTKFile("mesh1.vtk")
f2 = ot.SymbolicFunction(["x", "y", "z"], ["z^4-y^4+2*(y^2+z^2)-16"])
level2 = ot.LevelSet(f2, ot.LessOrEqual(), 0.0)
mesh2 = ot.LevelSetMesher([N] * 3).build(level2, ot.Interval([-5.0] * 3, [5.0] * 3))
# mesh2.exportToVTKFile("mesh2.vtk")
inter12 = otmeshing.IntersectionMesher().build([mesh1, mesh2])
volume = inter12.getVolume()
print("inter volume=", volume)
ott.assert_almost_equal(volume, 347.2445)
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
