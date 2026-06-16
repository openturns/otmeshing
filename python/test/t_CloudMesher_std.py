#! /usr/bin/env python

import math
import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

# basic 2d triangulation - BASIC method
mesher = otmeshing.CloudMesher()
print("mesher=", mesher)
assert "CloudMesher" in repr(mesher)
assert mesher.getClassName() == "CloudMesher"

vertices = [[3.0, 0.0], [2.0, 0.0], [2.0, 0.75], [2.5, 0.75], [3.0, 0.2]]
triangulation = mesher.build(vertices)
vol = triangulation.getVolume()
print(f"-- 2d BASIC={repr(triangulation)} vol={vol}")
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == 3
assert triangulation.isValid()
ott.assert_almost_equal(vol, 0.6125)
assert triangulation.isConvex()
# check all vertices are part of the input
verts = [[v[0], v[1]] for v in triangulation.getVertices()]
for v in verts:
    assert v in vertices, f"vertex {v} not in input"

# basic 2d triangulation - DELAUNAY method
mesher = otmeshing.CloudMesher(otmeshing.CloudMesher.DELAUNAY)
triangulation = mesher.build(vertices)
vol = triangulation.getVolume()
print(f"-- 2d DELAUNAY={repr(triangulation)} vol={vol}")
assert triangulation.getDimension() == 2
assert len(triangulation.getSimplices()) == 3
assert triangulation.isValid()
ott.assert_almost_equal(vol, 0.6125)
assert triangulation.isConvex()

# 1D case
mesher = otmeshing.CloudMesher()
triangulation = mesher.build([[2.5], [1.5], [3.0]])
assert triangulation.getDimension() == 1
assert triangulation.isValid()
assert triangulation.isConvex()
ott.assert_almost_equal(triangulation.getVolume(), 1.5)
# 1D special case: only min and max become vertices
assert len(triangulation.getVertices()) == 2
assert len(triangulation.getSimplices()) == 1
verts = triangulation.getVertices()
ott.assert_almost_equal(verts[0], [1.5])
ott.assert_almost_equal(verts[1], [3.0])

# nd triangulation of the unit hypercube
for method in [otmeshing.CloudMesher.BASIC, otmeshing.CloudMesher.DELAUNAY]:
    mesher = otmeshing.CloudMesher(method)
    for dim in range(1, 7):
        print(f"-- cube dim={dim} method={method}")
        vertices = ot.Box([0] * dim).generate()
        triangulation = mesher.build(vertices)
        assert triangulation.getDimension() == dim
        vol = triangulation.getVolume()
        print(f"vol={vol} triangulation={repr(triangulation)[:2000]}")
        assert len(triangulation.getVertices()) == len(vertices)
        assert triangulation.isValid()
        ott.assert_almost_equal(vol, 1.0)
        assert triangulation.isConvex()

# nd triangulation of the unit hypersphere
for method in [otmeshing.CloudMesher.BASIC, otmeshing.CloudMesher.DELAUNAY]:
    mesher = otmeshing.CloudMesher(method)
    for dim in range(1, 5):
        print(f"-- sphere dim={dim} method={method}")
        vertices = ot.Normal(dim).getSample(1000)
        for i in range(len(vertices)):
            vI = vertices[i]
            vertices[i] = vI / vI.norm()
        triangulation = mesher.build(vertices)
        assert triangulation.getDimension() == dim
        vol = triangulation.getVolume()
        print(f"vol={vol} triangulation={repr(triangulation)[:2000]}")
        if dim > 1:
            assert len(triangulation.getVertices()) == len(vertices)
        assert triangulation.isValid()
        vol_ref = math.pi ** (dim / 2) / math.gamma(dim / 2 + 1)
        ott.assert_almost_equal(vol, vol_ref, 0.1, 0.0)
        assert triangulation.isConvex()

# error cases: null dimension
with ott.assert_raises(TypeError):
    otmeshing.CloudMesher().build([])

# error cases: invalid empty point
with ott.assert_raises(TypeError):
    otmeshing.CloudMesher().build([[]])

# error cases: insufficient points for dimension
with ott.assert_raises(TypeError):
    otmeshing.CloudMesher().build([[1.0]])
with ott.assert_raises(TypeError):
    otmeshing.CloudMesher().build([[1.0, 2.0]])
