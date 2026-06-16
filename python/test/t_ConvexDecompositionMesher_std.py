#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import os
import otmeshing

ot.TESTPREAMBLE()

# Default constructor and __repr__
mesher = otmeshing.ConvexDecompositionMesher()
print("mesher=", mesher)
assert mesher.__repr__() == "class=ConvexDecompositionMesher"
assert not mesher.getUseSimplicesDecomposition()

# Copy constructor
mesher_copy = otmeshing.ConvexDecompositionMesher(mesher)
assert mesher_copy.__repr__() == mesher.__repr__()
assert (
    mesher_copy.getUseSimplicesDecomposition() == mesher.getUseSimplicesDecomposition()
)

# Save/load
study_file = "cdm_test.xml"
study = ot.Study(study_file)
study.add("mesher", mesher)
study.save()
loaded = otmeshing.ConvexDecompositionMesher()
study2 = ot.Study(study_file)
study2.load()
study2.fillObject("mesher", loaded)
assert loaded.getUseSimplicesDecomposition() == mesher.getUseSimplicesDecomposition()
os.remove(study_file)

# 2d snake + box
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
mesh1 = otmeshing.PolygonMesher().build(polyline)
mesh2 = ot.IntervalMesher([1] * 2).build(ot.Interval([-2] * 2, [-1] * 2))
mesh = otmeshing.UnionMesher().build([mesh1, mesh2])
print(repr(mesh))
decomposition = mesher.build(mesh)
assert len(decomposition) == 7
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    print(i, repr(convex), convex.getVolume())
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, mesh.getVolume())

# 3d surface mesh
vertices = [
    [0, 0, 0],
    [2, 0, 0],
    [2, 2, 0],
    [0, 2, 0],  # 0-3 bottom
    [0, 0, 2],
    [2, 0, 2],
    [0, 2, 2],  # 4-7 top
    [1, 1, 1],
    [2, 1, 1],
    [2, 2, 1],
    [1, 2, 1],  # notch bottom
    [1, 1, 2],
    [2, 1, 2],
    [1, 2, 2],
]  # notch top
simplices = [
    [0, 1, 5, 5],
    [0, 5, 4, 4],  # y=0
    [0, 4, 6, 6],
    [0, 6, 3, 3],  # x=0
    [0, 2, 1, 1],
    [0, 3, 2, 2],  # z=0
    [10, 9, 2, 2],
    [10, 2, 3, 3],
    [10, 3, 6, 6],
    [10, 6, 13, 13],  # y=2 wall
    [8, 2, 9, 9],
    [8, 1, 2, 2],
    [8, 5, 1, 1],
    [8, 12, 5, 5],  # x=2 wall
    [11, 13, 6, 6],
    [11, 6, 4, 4],
    [11, 4, 5, 5],
    [11, 5, 12, 12],  # z=2 wall
    [8, 7, 12, 12],
    [7, 11, 12, 12],  # y=1 notch
    [7, 10, 13, 13],
    [7, 13, 11, 11],  # x=1 notch
    [7, 8, 9, 9],
    [7, 9, 10, 10],
]  # z=1 notch
mesh = ot.Mesh(vertices, simplices)
print(mesh)
print(mesh.getVolume())
assert mesh.isValid()
assert not otmeshing.ConvexDecompositionMesher.IsConvex(mesh)

# build decomposition
decomposition = mesher.build(mesh)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    print(i, repr(convex), convex.getVolume())
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, 7.0)

# 3d volumetric mesh (two overlapping cubes)
vertices = (
    ot.IntervalMesher([1, 1, 1]).build(ot.Interval([0.0] * 3, [2.0] * 3)).getVertices()
)
vertices.add(
    ot.IntervalMesher([1, 1, 1]).build(ot.Interval([1.0] * 3, [3.0] * 3)).getVertices()
)
simplices = [
    [0, 1, 5, 7],
    [0, 3, 1, 7],
    [0, 5, 4, 7],
    [0, 4, 6, 7],
    [0, 6, 2, 7],
    [0, 2, 3, 7],
    [8, 9, 13, 15],
    [8, 11, 9, 15],
    [8, 13, 12, 15],
    [8, 12, 14, 15],
    [8, 14, 10, 15],
    [8, 10, 11, 15],
]
mesh_3d_overlap = ot.Mesh(vertices, simplices)
print(mesh_3d_overlap)
print(mesh_3d_overlap.getVolume())
assert mesh_3d_overlap.isValid()

# build decomposition
decomposition = mesher.build(mesh_3d_overlap)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    print(i, repr(convex), convex.getVolume())
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, 15.0)

# 3d disconnected cubes
mesh1 = ot.IntervalMesher([1] * 3).build(ot.Interval([-2.0] * 3, [-1.0] * 3))
mesh2 = ot.IntervalMesher([1] * 3).build(ot.Interval([1.0] * 3, [2.0] * 3))
mesh = otmeshing.UnionMesher().build([mesh1, mesh2])
decomposition = mesher.build(mesh)
assert len(decomposition) == 2
for convex in decomposition:
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)

# Create a 4-D torus
f = ot.SymbolicFunction(
    ["x0", "x1", "x2", "x3"], ["(x0^2 + x1^2 + x2^2 + x3^2 + 3)^2 - 16 * (x0^2 + x1^2)"]
)
levelSet = ot.LevelSet(f, ot.LessOrEqual(), 0.0)
N = 11
mesh_4d = ot.LevelSetMesher([N] * 4).build(
    levelSet, ot.Interval([-3.0] * 2 + [-1.0] * 2, [3.0] * 2 + [1.0] * 2)
)
print(mesh_4d)
print(mesh_4d.getVolume())
assert mesh_4d.isValid()

# build decomposition
decomposition = mesher.build(mesh_4d)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, mesh_4d.getVolume())

# Create a 2-D torus, i.e two disks
f = ot.SymbolicFunction(["x0", "x1"], ["(x0^2 + x1^2 + 3)^2 - 16 * (x0^2)"])
levelSet = ot.LevelSet(f, ot.LessOrEqual(), 0.0)
N = 21
mesh = ot.LevelSetMesher([N] * 2).build(
    levelSet, ot.Interval([-3.0] + [-1.0], [3.0] + [1.0])
)
print(mesh)
print(mesh.getVolume())
assert mesh.isValid()

# build decomposition
decomposition = mesher.build(mesh)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, mesh.getVolume())

# useSimplicesDecomposition flag in 2D
mesher_simple = otmeshing.ConvexDecompositionMesher()
mesher_simple.setUseSimplicesDecomposition(True)
decomposition = mesher_simple.build(mesh)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    volume_sum += convex.getVolume()
    assert convex.getSimplicesNumber() == 1
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, mesh.getVolume())

# Already-convex mesh in 2D returns a single component
convexMesh = ot.IntervalMesher([1, 1]).build(ot.Interval([0.0] * 2, [1.0] * 2))
assert otmeshing.ConvexDecompositionMesher.IsConvex(convexMesh)
decomposition = mesher.build(convexMesh)
assert len(decomposition) == 1, "2D convex mesh should decompose to 1 component"
assert otmeshing.ConvexDecompositionMesher.IsConvex(decomposition[0])

# Already-convex mesh in 3D returns a single component
convexMesh3D = ot.IntervalMesher([1, 1, 1]).build(ot.Interval([0.0] * 3, [1.0] * 3))
assert otmeshing.ConvexDecompositionMesher.IsConvex(convexMesh3D)
decomposition = mesher.build(convexMesh3D)
assert len(decomposition) == 1, "3D convex mesh should decompose to 1 component"
assert otmeshing.ConvexDecompositionMesher.IsConvex(decomposition[0])

# 3D volumetric mesh with useSimplicesDecomposition=True
decomposition = mesher_simple.build(mesh_3d_overlap)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, mesh_3d_overlap.getVolume())

# 4D with useSimplicesDecomposition=True
decomposition = mesher_simple.build(mesh_4d)
volume_sum = 0.0
for i, convex in enumerate(decomposition):
    volume_sum += convex.getVolume()
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
ott.assert_almost_equal(volume_sum, mesh_4d.getVolume())

# Flag get/set roundtrip
mesher3 = otmeshing.ConvexDecompositionMesher()
assert not mesher3.getUseSimplicesDecomposition()
mesher3.setUseSimplicesDecomposition(True)
assert mesher3.getUseSimplicesDecomposition()
mesher3.setUseSimplicesDecomposition(False)
assert not mesher3.getUseSimplicesDecomposition()

# 2D disconnected components
print("--- 2D disconnected components ---")
mesh1_2d = ot.IntervalMesher([1, 1]).build(ot.Interval([-2.0] * 2, [-1.0] * 2))
mesh2_2d = ot.IntervalMesher([1, 1]).build(ot.Interval([1.0] * 2, [2.0] * 2))
mesh_2d_disc = otmeshing.UnionMesher().build([mesh1_2d, mesh2_2d])
decomposition = mesher.build(mesh_2d_disc)
volume_sum = sum(c.getVolume() for c in decomposition)
ott.assert_almost_equal(volume_sum, mesh_2d_disc.getVolume())
assert len(decomposition) == 2
for convex in decomposition:
    assert otmeshing.ConvexDecompositionMesher.IsConvex(convex)
print("OK")

# Error cases: empty mesh should return empty or raise
print("--- Empty mesh (single vertex) ---")
try:
    empty = ot.Mesh([[0.0]])
    result = mesher.build(empty)
    assert len(result) == 0
    print("OK: empty result")
except Exception as e:
    print("OK: raised", type(e).__name__)
