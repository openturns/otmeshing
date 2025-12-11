#! /usr/bin/env python

import openturns as ot
import otmeshing

ot.TESTPREAMBLE()

# 2d triangulation
vertices = [[3.0, 0.0], [2.0, 0.0], [2.0, 0.75,], [2.5, 0.75], [3.0, 0.2]]
mesher = otmeshing.CloudMesher()
print("mesher=", mesher)
triangulation = mesher.build(vertices)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 2
assert triangulation.getVertices() == vertices
assert len(triangulation.getSimplices()) == 3
assert triangulation.getSimplices()[0] == [1, 0, 4]
assert triangulation.getSimplices()[1] == [1, 3, 2]
assert triangulation.getSimplices()[2] == [1, 4, 3]
assert triangulation.isValid()

# 3d triangulation of the unit cube
vertices = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
            [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
triangulation = mesher.build(vertices)
print("triangulation=", repr(triangulation))
assert triangulation.getDimension() == 3
assert triangulation.getVertices() == vertices
assert len(triangulation.getSimplices()) == 6
assert triangulation.isValid()

# nd triangulation of a random sample
for d in range(1, 8):
    vertices = ot.Normal(d).getSample(20)
    triangulation = mesher.build(vertices)
    print(f"d={d} triangulation={repr(triangulation)[:2000]}")
    assert triangulation.getDimension() == d
    assert triangulation.isValid()
