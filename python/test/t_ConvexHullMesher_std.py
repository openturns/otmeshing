#! /usr/bin/env python

import math
import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.ConvexHullMesher()
print("mesher=", mesher)

# check string representation
assert str(mesher) == "class=ConvexHullMesher"

# 1. empty sample -> error
with ott.assert_raises(TypeError):
    mesher.build(ot.Sample(0, 0))

# 2. not enough points
with ott.assert_raises(TypeError):
    mesher.build(ot.Sample([[0.0]]))
with ott.assert_raises(TypeError):
    mesher.build(ot.Sample([[0.0, 0.0], [1.0, 0.0]]))

# 3. 1D: segment [0.5, 1.5]
p = ot.Sample([[0.5], [1.5]])
hull = mesher.build(p)
assert hull.getDimension() == 1
assert hull.getIntrinsicDimension() == 1
assert hull.getVerticesNumber() == 2
assert hull.getSimplicesNumber() == 1
ott.assert_almost_equal(hull.getVolume(), 1.0)
assert hull.isValid()
assert hull.isConvex()

# 4. 1D: many collinear points -> hull is [min, max]
p = ot.Sample([[0.5], [1.5], [2.0], [3.0], [6.0]])
hull = mesher.build(p)
assert hull.getVerticesNumber() == 2
assert hull.getSimplicesNumber() == 1
ott.assert_almost_equal(hull.getVolume(), 5.5)
assert hull.isValid()
assert hull.isConvex()

# 5. 2D: triangle
p = ot.Sample([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
hull = mesher.build(p)
assert hull.getDimension() == 2
assert hull.getIntrinsicDimension() == 1
assert hull.getVerticesNumber() == 3
assert hull.getSimplicesNumber() == 3
ott.assert_almost_equal(hull.getVolume(), 2.0 + math.sqrt(2.0))
assert hull.isValid()
assert hull.isConvex()

# 6. 2D: square -- all 4 vertices on hull
p = ot.Sample([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
hull = mesher.build(p)
assert hull.getDimension() == 2
assert hull.getIntrinsicDimension() == 1
assert hull.getVerticesNumber() == 4
assert hull.getSimplicesNumber() == 4
ott.assert_almost_equal(hull.getVolume(), 4.0)
assert hull.isValid()
assert hull.isConvex()

# 7. 2D: square with interior points -> hull unchanged
p = ot.Sample(
    [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5], [0.2, 0.3], [0.8, 0.1]]
)
hull = mesher.build(p)
assert hull.getVerticesNumber() == 4
assert hull.getSimplicesNumber() == 4
ott.assert_almost_equal(hull.getVolume(), 4.0)
assert hull.isValid()
assert hull.isConvex()

# 8. 3D: tetrahedron -- all 4 vertices on hull
p = ot.Sample([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
hull = mesher.build(p)
assert hull.getDimension() == 3
assert hull.getIntrinsicDimension() == 2
assert hull.getVerticesNumber() == 4
assert hull.getSimplicesNumber() == 4
# 4 triangular faces: areas are 0.5, 0.5, 0.5, sqrt(3)/2 ~= 0.8660
ott.assert_almost_equal(hull.getVolume(), 1.5 + 0.5 * math.sqrt(3.0))
assert hull.isValid()
assert hull.isConvex()

# 9. 3D: cube -- critical test for Bug 1 (triggers "Qt" centroid in Qhull)
p = ot.Sample(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
)
hull = mesher.build(p)
assert hull.getDimension() == 3
assert hull.getIntrinsicDimension() == 2
# CGAL: 8 vertices; Qhull: 14 (8 corners + 6 face centroids from "Qt")
assert hull.getVerticesNumber() >= 8
# 6 faces x 2 triangles = 12 simplices (both backends agree)
assert hull.getSimplicesNumber() == 12
ott.assert_almost_equal(hull.getVolume(), 6.0)
assert hull.isValid()
assert hull.isConvex()

# 10. 3D: cube with interior point -> hull vertices unchanged
p.add(ot.Point([0.5, 0.5, 0.5]))
hull = mesher.build(p)
assert hull.getVerticesNumber() >= 8
assert hull.getSimplicesNumber() == 12
ott.assert_almost_equal(hull.getVolume(), 6.0)

# 11. 4D: simplex (5 vertices)
p = ot.Sample(
    [
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
hull = mesher.build(p)
assert hull.getDimension() == 4
assert hull.getIntrinsicDimension() == 3
assert hull.getVerticesNumber() == 5
assert hull.getSimplicesNumber() == 5
ott.assert_almost_equal(hull.getVolume(), 1.0)
assert hull.isValid()
assert hull.isConvex()

# 12. 4D: hypercube (16 vertices)
p = ot.Box([0] * 4).generate()
hull = mesher.build(p)
assert hull.getDimension() == 4
assert hull.getIntrinsicDimension() == 3
assert hull.getVerticesNumber() >= 16
ott.assert_almost_equal(hull.getVolume(), 8.0)
assert hull.isValid()
assert hull.isConvex()

# 13. reusing the same mesher for multiple builds
hull1 = mesher.build(ot.Sample([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]))
hull2 = mesher.build(ot.Sample([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]))
assert hull1.getVerticesNumber() == 3
assert hull2.getVerticesNumber() == 3

# 14. ND: Gaussian random sample (dim 1-4)
for dim in range(1, 5):
    print(f"-- gaussian dim={dim}")
    p = ot.Normal(dim).getSample(1000)
    hull = mesher.build(p)
    assert hull.getDimension() == dim
    if dim > 1:
        assert hull.getIntrinsicDimension() == dim - 1
        assert hull.isValid()
    # Gaussian cloud has many interior points; hull has fewer
    assert hull.getVerticesNumber() < p.getSize()
