#! /usr/bin/env python

import math
import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

# basic 2d triangulation
vertices = [[3.0, 0.0], [2.0, 0.0], [2.0, 0.75,], [2.5, 0.75], [3.0, 0.2]]
mesher = otmeshing.CloudMesher()
print("mesher=", mesher)
triangulation = mesher.build(vertices)
vol = triangulation.getVolume()
print(f"2d triangulation={repr(triangulation)} vol={vol}")
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == 3
assert triangulation.isValid()
ott.assert_almost_equal(vol, 0.6125)

# nd triangulation of the unit hypercube
for dim in range(1, 9):
    vertices = ot.Box([0] * dim).generate()
    triangulation = mesher.build(vertices)
    vol = triangulation.getVolume()
    print(f"cube dim={dim} vol={vol} triangulation={repr(triangulation)[:2000]}", flush=True)
    assert triangulation.getDimension() == dim
    if dim > 1:
        assert len(triangulation.getVertices()) == len(vertices)
    assert triangulation.isValid()
    ott.assert_almost_equal(vol, 1.0)

# nd triangulation of the unit hypersphere
for dim in range(1, 5):
    vertices = ot.Normal(dim).getSample(1000)
    for i in range(len(vertices)):
        vI = vertices[i]
        vertices[i] = vI / vI.norm()
    triangulation = mesher.build(vertices)
    vol = triangulation.getVolume()
    print(f"sphere dim={dim} vol={vol} triangulation={repr(triangulation)[:2000]}")
    assert triangulation.getDimension() == dim
    if dim > 1:
        assert len(triangulation.getVertices()) == len(vertices)
    assert triangulation.isValid()
    vol_ref = math.pi**(dim / 2) / math.gamma(dim / 2 + 1)
    ott.assert_almost_equal(vol, vol_ref, 0.1, 0.0)
