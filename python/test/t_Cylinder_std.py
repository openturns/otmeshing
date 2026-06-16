#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing as otm
import os

ot.TESTPREAMBLE()

# Test default constructor
print("--- Default constructor ---")
cyl_default = otm.Cylinder()
assert cyl_default.getDimension() == 0
print("OK")

# Test parameter constructor with a simple square base (exact geometry)
print("--- Cylinder with square base ---")
base_mesh = ot.IntervalMesher([5, 5]).build(ot.Interval([-1.0, -1.0], [1.0, 1.0]))
extension = ot.Interval([-2.0], [2.0])
injection = [2]
M = 3
cyl_square = otm.Cylinder(base_mesh, extension, injection, M)

assert cyl_square.getDimension() == 3
ott.assert_almost_equal(
    cyl_square.getVolume(), base_mesh.getVolume() * extension.getVolume()
)

# check base
retrieved_base = cyl_square.getBase()
ott.assert_almost_equal(retrieved_base.getVertices(), base_mesh.getVertices())
assert retrieved_base.getSimplices() == base_mesh.getSimplices()

# check extension
retrieved_ext = cyl_square.getExtension()
assert retrieved_ext.getLowerBound()[0] == -2.0
assert retrieved_ext.getUpperBound()[0] == 2.0

# check injection
retrieved_inj = cyl_square.getInjection()
assert retrieved_inj == injection

# check discretization
assert cyl_square.getDiscretization() == M
print("OK")

# Test getVertices
print("--- getVertices ---")
ext_vertices = ot.IntervalMesher([M]).build(extension).getVertices()
vertices = cyl_square.getVertices()
assert vertices.getSize() == base_mesh.getVerticesNumber() * ext_vertices.getSize()
assert vertices.getDimension() == 3
print("OK")

# Test getBoundingBox
print("--- getBoundingBox ---")
bbox = cyl_square.getBoundingBox()
ott.assert_almost_equal(bbox.getLowerBound(), [-1.0, -1.0, -2.0])
ott.assert_almost_equal(bbox.getUpperBound(), [1.0, 1.0, 2.0])
print("OK")

# Test isConvex
print("--- isConvex ---")
assert cyl_square.isConvex()
print("OK")

# Test computeMesh
print("--- computeMesh ---")
computed_mesh = cyl_square.computeMesh()
assert computed_mesh.getVerticesNumber() > 0
print("OK")

# Test string representation
print("--- __repr__ / __str__ ---")
assert cyl_square.__repr__() == "class=Cylinder"
assert cyl_square.__str__() == "class=Cylinder"
print("OK")


# Test save/load
print("--- save/load ---")

study_file = "cylinder_test.xml"
study = ot.Study(study_file)
study.add("cylinder", cyl_square)
study.save()

loaded = otm.Cylinder()
study2 = ot.Study(study_file)
study2.load()
study2.fillObject("cylinder", loaded)
assert loaded.getDimension() == cyl_square.getDimension()
ott.assert_almost_equal(loaded.getVolume(), cyl_square.getVolume())
ott.assert_almost_equal(
    loaded.getBoundingBox().getLowerBound(), cyl_square.getBoundingBox().getLowerBound()
)
ott.assert_almost_equal(
    loaded.getBoundingBox().getUpperBound(), cyl_square.getBoundingBox().getUpperBound()
)
assert loaded.getDiscretization() == cyl_square.getDiscretization()
os.remove(study_file)
print("OK")

# Test non-convex cylinder
print("--- Non-convex cylinder ---")
snake_polyline = [
    [0.0, 0.0],
    [0.0, 5.0],
    [6.0, 5.0],
    [6.0, 0.0],
    [2.0, 0.0],
    [2.0, 3.0],
    [4.0, 3.0],
    [4.0, 2.0],
    [3.0, 2.0],
    [3.0, 1.0],
    [5.0, 1.0],
    [5.0, 4.0],
    [1.0, 4.0],
    [1.0, 0.0],
]
snake_mesh = otm.PolygonMesher().build(snake_polyline)
nonconvex_cyl = otm.Cylinder(snake_mesh, extension, injection, M)
assert nonconvex_cyl.getDimension() == 3
assert not nonconvex_cyl.isConvex()
ott.assert_almost_equal(
    nonconvex_cyl.getVolume(), snake_mesh.getVolume() * extension.getVolume()
)
print("OK")

# Test 1D base (extension in 2D)
print("--- 1D base cylinder ---")
base_1d = ot.IntervalMesher([5]).build(ot.Interval([-1.0], [1.0]))
ext_2d = ot.Interval([-1.0] * 2, [1.0] * 2)
inj_2d = [1, 2]
M2 = 2
cyl_1d = otm.Cylinder(base_1d, ext_2d, inj_2d, M2)
assert cyl_1d.getDimension() == 3
ott.assert_almost_equal(cyl_1d.getVolume(), base_1d.getVolume() * ext_2d.getVolume())
assert cyl_1d.isConvex()
bbox_1d = cyl_1d.getBoundingBox()
ott.assert_almost_equal(bbox_1d.getLowerBound(), [-1.0, -1.0, -1.0])
ott.assert_almost_equal(bbox_1d.getUpperBound(), [1.0, 1.0, 1.0])
print("OK")

# Test error cases
print("--- Error cases ---")
# Invalid injection size
with ott.assert_raises(TypeError):
    otm.Cylinder(base_mesh, extension, [2, 3], M)
print("Caught invalid injection size OK")

# Injection index out of range
with ott.assert_raises(TypeError):
    otm.Cylinder(base_mesh, extension, [5], M)
print("Caught invalid injection index OK")

# Keep existing disc-based tests
print("--- Disc-based cylinder ---")
N = 20
dim = 3
xc1 = 0.0
yc1 = 0.0
R1 = 1.0
H1 = 4.0
f1 = ot.SymbolicFunction(["x", "y"], [f"(x-({xc1}))^2+(y-({yc1}))^2"])
levelSet1 = ot.LevelSet(f1, ot.LessOrEqual(), R1**2)
base1 = ot.LevelSetMesher([N] * 2).build(
    levelSet1, ot.Interval([xc1 - R1, yc1 - R1], [xc1 + R1, yc1 + R1])
)
extension1 = ot.Interval([-H1 / 2] * (dim - 2), [H1 / 2] * (dim - 2))
injection1 = list(range(2, dim))
cyl1 = otm.Cylinder(base1, extension1, injection1, M)
assert cyl1.getDimension() == dim
assert cyl1.isConvex()
# volume = base_area * 4 (extension length)
vol1 = cyl1.getVolume()
ott.assert_almost_equal(vol1, base1.getVolume() * extension1.getVolume())
print("Cylinder 1 volume:", vol1)

print("--- Disc-based cylinder 2 ---")
yc2 = 0.8
zc2 = 0.0
R2 = 0.5
H2 = 4.0
f2 = ot.SymbolicFunction(["y", "z"], [f"(y-({yc2}))^2+(z-({zc2}))^2"])
levelSet2 = ot.LevelSet(f2, ot.LessOrEqual(), R2**2)
base2 = ot.LevelSetMesher([N] * 2).build(
    levelSet2, ot.Interval([yc2 - R2, zc2 - R2], [yc2 + R2, zc2 + R2])
)
extension2 = ot.Interval([-H2 / 2] * (dim - 2), [H2 / 2] * (dim - 2))
injection2 = [0] + list(range(3, dim))
cyl2 = otm.Cylinder(base2, extension2, injection2, M)
assert cyl2.isConvex()
vol2 = cyl2.getVolume()
ott.assert_almost_equal(vol2, base2.getVolume() * extension2.getVolume())
print("Cylinder 2 volume:", vol2)

# Test computeMesh on disc cylinders
print("--- computeMesh on disc cylinders ---")
mesh1 = cyl1.computeMesh()
assert mesh1.getVerticesNumber() > 0
assert mesh1.getDimension() == dim
assert mesh1.isValid()
print("Mesh 1 vertices:", mesh1.getVerticesNumber())

mesh2 = cyl2.computeMesh()
assert mesh2.getVerticesNumber() > 0
assert mesh2.getDimension() == dim
assert mesh2.isValid()
print("Mesh 2 vertices:", mesh2.getVerticesNumber())

# Intersect cylinders (existing test)
print("--- Intersect cylinders ---")
inter12 = otm.IntersectionMesher().buildConvex([mesh1, mesh2])
print("Intersection volume:", inter12.getVolume())
ott.assert_almost_equal(inter12.getVolume(), 0.778142671)

print("All tests passed")
