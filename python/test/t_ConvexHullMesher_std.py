#! /usr/bin/env python

import openturns as ot
import otmeshing

ot.TESTPREAMBLE()

# basic 2d convex hull
mesher = otmeshing.ConvexHullMesher()
print("mesher=", mesher)
vertices = [[0.0, 0.0], [0.0, 3.0], [1.0, 2.0], [0.5, 0.5], [1.0, 0.5], [2.0, 2.0], [2.0, 0.0]]
hull = mesher.build(vertices)
print(f"-- 2d hull={repr(hull)}")
assert hull.getDimension() == 2
assert hull.getIntrinsicDimension() == 1
assert len(hull.getVertices()) == 4
assert len(hull.getSimplices()) == 4
assert hull.isValid()

# nd hull of the unit hypercube
for dim in range(1, 7):
    print(f"-- cube dim={dim}")
    vertices = ot.Box([0] * dim).generate()
    hull = mesher.build(vertices)
    print(hull.getDimension(), hull.getIntrinsicDimension())
    print(f"hull={repr(hull)[:2000]}")
    assert hull.getDimension() == dim
    if dim > 1:
        assert hull.getIntrinsicDimension() == dim - 1
        assert len(hull.getVertices()) == len(vertices)
    assert hull.isValid()

# nd hull of a Gaussian sample
for dim in range(1, 5):
    print(f"-- sphere dim={dim}")
    vertices = ot.Normal(dim).getSample(1000)
    hull = mesher.build(vertices)
    print(f"hull={repr(hull)[:2000]}")
    assert hull.getDimension() == dim
    if dim > 1:
        assert hull.getIntrinsicDimension() == dim - 1
        assert len(hull.getVertices()) < len(vertices)
    assert hull.isValid()
