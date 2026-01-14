#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

# 3d surface mesh
vertices = [
    [0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0],  # 0–3 bottom
    [0, 0, 2], [2, 0, 2], [0, 2, 2],  # 4–7 top
    [1, 1, 1], [2, 1, 1], [2, 2, 1], [1, 2, 1],  # notch bottom
    [1, 1, 2], [2, 1, 2], [1, 2, 2]]  # notch top
simplices = [
    [0, 1, 5, 5], [0, 5, 4, 4],  # y=0
    [0, 4, 6, 6], [0, 6, 3, 3],  # x=0
    [0, 2, 1, 1], [0, 3, 2, 2],  # z=0
    [10, 9, 2, 2], [10, 2, 3, 3], [10, 3, 6, 6], [10, 6, 13, 13],  # y=2 wall
    [8, 2, 9, 9], [8, 1, 2, 2], [8, 5, 1, 1], [8, 12, 5, 5],  # x=2 wall
    [11, 13, 6, 6], [11, 6, 4, 4], [11, 4, 5, 5], [11, 5, 12, 12],  # z=2 wall
    [8, 7, 12, 12], [7, 11, 12, 12],  # y=1 notch
    [7, 10, 13, 13], [7, 13, 11, 11],  # x=1 notch
    [7, 8, 9, 9], [7, 9, 10, 10]]  # z=1 notch
mesh = ot.Mesh(vertices, simplices)
print(mesh)
print(mesh.getVolume())
assert mesh.isValid()

# build decomposition
mesher = otmeshing.ConvexDecompositionMesher()
print("mesher=", mesher)
decomposition = mesher.build(mesh)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    print(i, repr(convex), convex.getVolume())
    volume_sum += convex.getVolume()
ott.assert_almost_equal(volume_sum, 7.0)

# 3d volumetric mesh (two overlapping cubes)
vertices = ot.IntervalMesher([1, 1, 1]).build(ot.Interval([0.0] * 3, [2.0] * 3)).getVertices()
vertices.add(ot.IntervalMesher([1, 1, 1]).build(ot.Interval([1.0] * 3, [3.0] * 3)).getVertices())
simplices = [[0, 1, 5, 7], [0, 3, 1, 7], [0, 5, 4, 7], [0, 4, 6, 7], [0, 6, 2, 7], [0, 2, 3, 7],
             [8, 9, 13, 15], [8, 11, 9, 15], [8, 13, 12, 15], [8, 12, 14, 15], [8, 14, 10, 15], [8, 10, 11, 15]]
mesh = ot.Mesh(vertices, simplices)
print(mesh)
print(mesh.getVolume())
assert mesh.isValid()

# build decomposition
decomposition = mesher.build(mesh)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    print(i, repr(convex), convex.getVolume())
    volume_sum += convex.getVolume()
ott.assert_almost_equal(volume_sum, 15.0)
