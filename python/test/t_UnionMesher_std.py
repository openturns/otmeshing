#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing

ot.TESTPREAMBLE()

mesher = otmeshing.UnionMesher()
print(mesher)
print(repr(mesher))

# 1. Empty collection -> empty mesh
empty = mesher.build([])
print(f"Empty collection: {empty}")
ott.assert_almost_equal(empty.getDimension(), 0)
ott.assert_almost_equal(empty.getVerticesNumber(), 0)

# 2. Single mesh -> returns compressed mesh (no-op for valid input)
for dim in range(2, 5):
    mesh = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
    single = mesher.build([mesh])
    print(f"Single dim={dim}: {single}")
    assert single.isValid()
    ott.assert_almost_equal(single.getVerticesNumber(), mesh.getVerticesNumber())
    ott.assert_almost_equal(single.getSimplicesNumber(), mesh.getSimplicesNumber())
    ott.assert_almost_equal(single.getVolume(), 1.0)

# 3. Non-overlapping meshes (original test, extended)
for dim in range(2, 6):
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))
    union = mesher.build([mesh1, mesh2])
    print(f"Disjoint dim={dim}: {union}")
    assert union.isValid()
    ott.assert_almost_equal(union.getVolume(), 2.0)

# 4. Touching meshes (shared boundary) -- exercises vertex dedup
# mesh1 = [0,1]^dim, mesh2 = [1,2]x[0,1]^{dim-1}, sharing the face at x=1
for dim in range(2, 5):
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))  # [0,1]^dim
    lower = [1.0] + [0.0] * (dim - 1)
    upper = [2.0] + [1.0] * (dim - 1)
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval(lower, upper))
    union = mesher.build([mesh1, mesh2])
    print(f"Touching dim={dim}: {union}")
    assert union.isValid()
    ott.assert_almost_equal(union.getVolume(), 2.0)
    # shared (dim-1)-face has 2^{dim-1} vertices
    nVertices = union.getVerticesNumber()
    nNoDedup = 2 ** (dim + 1)
    nShared = 2 ** (dim - 1)
    expectedVertices = nNoDedup - nShared
    print(f"  vertices={nVertices} expected={expectedVertices}")
    ott.assert_almost_equal(nVertices, expectedVertices)

# 5. Three disjoint meshes combined
for dim in range(2, 4):
    mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
    mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))
    mesh3 = ot.IntervalMesher([1] * dim).build(ot.Interval([4.0] * dim, [5.0] * dim))
    union = mesher.build([mesh1, mesh2, mesh3])
    print(f"Three dim={dim}: {union}")
    assert union.isValid()
    ott.assert_almost_equal(union.getVolume(), 3.0)

# 6. CompressMesh static method: duplicate vertices
mesh = ot.IntervalMesher([1] * 2).build(ot.Interval(2))
compressed = otmeshing.UnionMesher.CompressMesh(mesh)
print(f"Already compressed: {compressed}")
assert compressed.isValid()
ott.assert_almost_equal(compressed.getVerticesNumber(), mesh.getVerticesNumber())
ott.assert_almost_equal(compressed.getSimplicesNumber(), mesh.getSimplicesNumber())
ott.assert_almost_equal(compressed.getVolume(), 1.0)

# 7. CompressMesh static method: collapse duplicates
vertices = ot.Sample(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ]
)
simplices = ot.IndicesCollection([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]])
mesh = ot.Mesh(vertices, simplices, False)
compressed = otmeshing.UnionMesher.CompressMesh(mesh)
print(f"Duplicates collapsed: {compressed}")
assert compressed.isValid()
ott.assert_almost_equal(compressed.getVerticesNumber(), 4)
ott.assert_almost_equal(compressed.getVolume(), 2.0)

# 8. CompressMesh static method: used-vertex-only output
# A mesh where only some vertices are referenced by simplices
vertices = ot.Sample(
    [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [42.0, 42.0]]
)  # unused vertex
simplices = ot.IndicesCollection([[0, 1, 2], [0, 2, 3]])
mesh = ot.Mesh(vertices, simplices, False)
compressed = otmeshing.UnionMesher.CompressMesh(mesh)
print(f"Unused dropped: {compressed}")
assert compressed.isValid()
ott.assert_almost_equal(compressed.getVerticesNumber(), 4)
ott.assert_almost_equal(compressed.getVolume(), 1.0)

# 9. build with a mesh that has an unused vertex
union = mesher.build([mesh])
print(f"Build with unused: {union}")
assert union.isValid()
ott.assert_almost_equal(union.getVerticesNumber(), 4)
ott.assert_almost_equal(union.getVolume(), 1.0)

# 10. Empty input mesh
empty_mesh = ot.Mesh(ot.Sample(0, 2), ot.IndicesCollection())
single = mesher.build([empty_mesh])
print(f"Empty mesh: {single}")
assert single.isValid()

# 11. Two empty meshes
union = mesher.build([empty_mesh, empty_mesh])
print(f"Two empty: {union}")
assert union.isValid()

# 12. Dimension mismatch should raise
mesh1 = ot.IntervalMesher([1] * 2).build(ot.Interval(2))
mesh2 = ot.IntervalMesher([1] * 3).build(ot.Interval(3))
with ott.assert_raises(TypeError):
    mesher.build([mesh1, mesh2])
print("Dimension mismatch correctly raised")

# 13. CompressMesh: duplicate vertices with non-zero range (tolerance > 0)
# Mix of unique and duplicate vertices so tolerance is non-zero
pts = [
    [0.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [0.0, 1.0],
    [0.0, 0.0],
    [1.0, 0.0],
]  # last two are duplicates of first two
simplices = ot.IndicesCollection([[0, 1, 2], [0, 2, 3], [4, 5, 2], [4, 2, 3]])
vertices = ot.Sample(pts)
mesh = ot.Mesh(vertices, simplices, False)
compressed = otmeshing.UnionMesher.CompressMesh(mesh)
print(f"Duplicates with non-zero range: {compressed}")
assert compressed.isValid()
# 4 unique vertices + 2 duplicates -> 4 compressed vertices
ott.assert_almost_equal(compressed.getVerticesNumber(), 4)
ott.assert_almost_equal(compressed.getVolume(), 2.0)
