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

# 14. Coincident meshes union to their common mesh (no double counting)
square = ot.IntervalMesher([1] * 2).build(ot.Interval(2))
coincident = mesher.build([square, square])
print(f"Coincident: {coincident}")
assert coincident.isValid()
ott.assert_almost_equal(coincident.getVolume(), 1.0)
ott.assert_almost_equal(coincident.getSimplicesNumber(), 2)
ott.assert_almost_equal(coincident.getVerticesNumber(), 4)


def edge_counts_2d(mesh):
    counts = {}
    for simplex in mesh.getSimplices():
        for e in range(3):
            key = tuple(sorted((simplex[e], simplex[(e + 1) % 3])))
            counts[key] = counts.get(key, 0) + 1
    return counts


# 15. T-junction: shared vertices connected differently across meshes
# A = [0,1]^2, B = [1,2]x[0,1] with a mid node (1, 0.5) on the shared edge
meshA = ot.Mesh([[0, 0], [1, 0], [1, 1], [0, 1]], [[0, 1, 2], [0, 2, 3]])
meshB = ot.Mesh(
    [[1, 0], [2, 0], [2, 1], [1, 1], [1, 0.5]], [[0, 1, 2], [0, 2, 4], [2, 3, 4]]
)
tjunction = mesher.build([meshA, meshB])
print(f"T-junction: {tjunction}")
assert tjunction.isValid()
ott.assert_almost_equal(tjunction.getVolume(), 2.0)
ott.assert_almost_equal(tjunction.getVerticesNumber(), 7)
ott.assert_almost_equal(tjunction.getSimplicesNumber(), 6)
# conforming interface: shared edge split node-for-node, each piece used twice
counts = edge_counts_2d(tjunction)
vertices = tjunction.getVertices()
interface = [
    k
    for k in counts
    if abs(vertices[k[0]][0] - 1.0) < 1e-12 and abs(vertices[k[1]][0] - 1.0) < 1e-12
]
assert len(interface) == 2, "shared edge must be split in two"
assert all(counts[k] == 2 for k in interface), "interface edges must be shared twice"

# 16. 3D T-junction: mismatched face resolutions across the shared face.
# A = [0,1]^3, B = [0,1]^2x[1,2] refined 2x2 on the shared face.
meshA3D = ot.IntervalMesher([1] * 3).build(ot.Interval([0.0] * 3, [1.0] * 3))
meshB3D = ot.IntervalMesher([2, 2, 1]).build(ot.Interval([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]))
tjunction3D = mesher.build([meshA3D, meshB3D])
print(f"T-junction 3D: {tjunction3D}")
assert tjunction3D.isValid()
ott.assert_almost_equal(tjunction3D.getVolume(), 2.0)
# no new vertices: splits reuse welded nodes only (8 + 18 - 4 corners)
ott.assert_almost_equal(tjunction3D.getVerticesNumber(), 22)
# refinement fired (6 + 24 base simplices subdivided)
assert tjunction3D.getSimplicesNumber() > 30

# 17. Reset-proof: module keys survive ResourceMap.Reset() (Sphinx plot
# pre-code calls Reset before every figure)
ot.ResourceMap.Reset()
resetUnion = mesher.build(
    [
        ot.IntervalMesher([1] * 2).build(ot.Interval(2)),
        ot.IntervalMesher([1] * 2).build(ot.Interval([2.0, 0.0], [3.0, 1.0])),
    ]
)
assert resetUnion.isValid()
ott.assert_almost_equal(resetUnion.getVolume(), 2.0)
