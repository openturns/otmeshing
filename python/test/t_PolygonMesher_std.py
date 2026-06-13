#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing
import math

ot.TESTPREAMBLE()

mesher = otmeshing.PolygonMesher()
print("mesher=", mesher)

# 2d triangulation of a convex polygon (regular n-sided)
polyline = []
n = 20
for i in range(n):
    r = 1.0
    theta = i * 2.0 * math.pi / n
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    polyline.append([x, y])
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == n - 2
assert triangulation.isValid()
# area of unit regular n-gon: n/2 * sin(2*pi/n)
area_ref = 0.5 * n * math.sin(2.0 * math.pi / n)
ott.assert_almost_equal(triangulation.getVolume(), area_ref)

# 2d triangulation of a non-convex polygon (snail-like)
polyline = [
    [0, 0],
    [0, 5],
    [6, 5],
    [6, 0],
    [2, 0],
    [2, 3],
    [4, 3],
    [4, 2],
    [3, 2],
    [3, 1],
    [5, 1],
    [5, 4],
    [1, 4],
    [1, 0],
]
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == 12
assert triangulation.isValid()
# area of this polygon: outer box 6x5 minus inner cutout 4x4 + some
# compute numerically: polygon area using shoelace formula
xs = [p[0] for p in polyline]
ys = [p[1] for p in polyline]
area_ref = 0.5 * abs(
    sum(
        xs[i] * ys[(i + 1) % len(polyline)] - xs[(i + 1) % len(polyline)] * ys[i]
        for i in range(len(polyline))
    )
)
ott.assert_almost_equal(triangulation.getVolume(), area_ref)

# simple triangle (minimal case)
polyline = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == 1
assert triangulation.isValid()
ott.assert_almost_equal(triangulation.getVolume(), 0.5)

# 2D polygon with points in clockwise order
polyline = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == 2
assert triangulation.isValid()
ott.assert_almost_equal(triangulation.getVolume(), 1.0)

# 3D coplanar points (z=0 plane)
polyline = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 3
assert len(triangulation.getSimplices()) == 2
assert triangulation.isValid()
ott.assert_almost_equal(triangulation.getVolume(), 0.0)

# check that the mesh simplices have the right size (dimension+1)
for simplex in triangulation.getSimplices():
    assert len(simplex) == 4  # 3D -> 4 vertices per simplex

# 5D coplanar points
polyline = [
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0, 0.0],
]
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 5
assert len(triangulation.getSimplices()) == 2
assert triangulation.isValid()

# error cases
# too few points
try:
    mesher.build([[0.0, 0.0], [1.0, 0.0]])
    assert False
except RuntimeError:
    pass

# redundant vertex
try:
    mesher.build([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 0.0]])
    assert False
except RuntimeError:
    pass

# collinear points (intrinsic dimension 1)
try:
    mesher.build([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    assert False
except RuntimeError:
    pass
