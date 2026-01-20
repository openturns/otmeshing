#! /usr/bin/env python

import openturns as ot
import otmeshing
import math

ot.TESTPREAMBLE()

# 2d triangulation of a convex polygon (regular n-sided)
polyline = []
n = 20
for i in range(n):
    r = 1.0
    theta = i * 2.0 * math.pi / n
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    polyline.append([x, y])
mesher = otmeshing.PolygonMesher()
print("mesher=", mesher)
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert triangulation.getVertices() == polyline
assert len(triangulation.getSimplices()) == n - 2
assert triangulation.isValid()

# 2d triangulation of a non-convex polygon (snail-like)
polyline = [[0, 0], [0, 5], [6, 5], [6, 0], [2, 0], [2, 3], [4, 3],
            [4, 2], [3, 2], [3, 1], [5, 1], [5, 4], [1, 4], [1, 0]]
triangulation = mesher.build(polyline)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert triangulation.getVertices() == polyline
assert len(triangulation.getSimplices()) == 12
assert triangulation.isValid()
