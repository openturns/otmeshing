#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import otmeshing as otm

ot.TESTPREAMBLE()

# 1. Default constructor
print("=" * 40)
print("Default constructor")
mesher = otm.FunctionGraphMesher()
print("mesher =", mesher)

# 2. 1D linear function -- exact volume with any discretization
print("=" * 40)
print("1D linear function")
f = ot.SymbolicFunction(["x"], ["x"])
f.setName("linear")
f.setInputDescription(["x"])
f.setOutputDescription(["y"])

inputInterval = ot.Interval([0.0], [1.0])
inputDiscretization = [5]
mesher = otm.FunctionGraphMesher(inputInterval, inputDiscretization)
print("mesher =", mesher)

# subGraph=True, outputIndex=1 (output last)
mesh = mesher.build(f, 1, -1.0, 2.0, 1, True)
print("mesh (subGraph=T, outIdx=1):", mesh)
print(
    "  dim=%d vert=%d simp=%d vol=%.12g"
    % (
        mesh.getDimension(),
        mesh.getVerticesNumber(),
        mesh.getSimplicesNumber(),
        mesh.getVolume(),
    )
)
assert mesh.getDimension() == 2
assert mesh.getVerticesNumber() == 12  # (5+1) * (1+1)
assert mesh.getSimplicesNumber() == 10  # 5*1*2
ott.assert_almost_equal(mesh.getVolume(), 1.5)  # int(x+1)dx = 3/2
assert mesh.getDescription() == ["x", "y"]

# subGraph=False, outputIndex=1
mesh2 = mesher.build(f, 1, -1.0, 2.0, 1, False)
print("mesh (subGraph=F, outIdx=1): vol=%.12g" % mesh2.getVolume())
ott.assert_almost_equal(mesh2.getVolume(), 1.5)  # int(2-x)dx = 3/2

# outputIndex=0 (output first)
mesh3 = mesher.build(f, 0, -1.0, 2.0, 1, True)
print(
    "mesh (subGraph=T, outIdx=0): dim=%d vol=%.12g"
    % (mesh3.getDimension(), mesh3.getVolume())
)
assert mesh3.getDimension() == 2
ott.assert_almost_equal(mesh3.getVolume(), 1.5)
assert mesh3.getDescription() == ["y", "x"]

# outputDiscretization > 1
mesh4 = mesher.build(f, 1, -1.0, 2.0, 3, True)
print(
    "mesh (subGraph=T, outDisc=3): vert=%d simp=%d vol=%.12g"
    % (mesh4.getVerticesNumber(), mesh4.getSimplicesNumber(), mesh4.getVolume())
)
assert mesh4.getVerticesNumber() == 24  # (5+1) * (3+1)
assert mesh4.getSimplicesNumber() == 30  # 5*3*2
ott.assert_almost_equal(mesh4.getVolume(), 1.5)

# 3. 1D constant function
print("=" * 40)
print("1D constant function")
f_zero = ot.SymbolicFunction(["x"], ["0.0"])
mesher_cst = otm.FunctionGraphMesher(ot.Interval([0.0], [1.0]), [1])
mesh = mesher_cst.build(f_zero, 1, -1.0, 1.0, 1, True)
print("  subgraph volume: %.12g" % mesh.getVolume())
ott.assert_almost_equal(mesh.getVolume(), 1.0)  # int(0-(-1))dx = 1
mesh = mesher_cst.build(f_zero, 1, -1.0, 1.0, 1, False)
print("  supergraph volume: %.12g" % mesh.getVolume())
ott.assert_almost_equal(mesh.getVolume(), 1.0)  # int(1-0)dx = 1

# 4. 2D linear function -- exact volume
print("=" * 40)
print("2D linear function")
f2d = ot.SymbolicFunction(["x0", "x1"], ["x0 + x1"])
f2d.setName("linear2d")
f2d.setInputDescription(["x", "y"])
f2d.setOutputDescription(["z"])

mesher3 = otm.FunctionGraphMesher(ot.Interval([0.0, 0.0], [1.0, 1.0]), [1, 1])
mesh = mesher3.build(f2d, 2, -1.0, 3.0, 1, True)
print(
    "mesh (subGraph=T, outIdx=2): dim=%d vert=%d simp=%d vol=%.12g"
    % (
        mesh.getDimension(),
        mesh.getVerticesNumber(),
        mesh.getSimplicesNumber(),
        mesh.getVolume(),
    )
)
assert mesh.getDimension() == 3
ott.assert_almost_equal(mesh.getVolume(), 2.0)  # int(x0+x1+1) = 2

# subGraph=False
mesh = mesher3.build(f2d, 2, -1.0, 3.0, 1, False)
print("  supergraph volume: %.12g" % mesh.getVolume())
ott.assert_almost_equal(mesh.getVolume(), 2.0)  # int(3-x0-x1) = 2

# outputIndex at alternative positions
for idx in [0, 1]:
    mesh = mesher3.build(f2d, idx, -1.0, 3.0, 1, True)
    print("  outputIndex=%d volume: %.12g" % (idx, mesh.getVolume()))
    ott.assert_almost_equal(mesh.getVolume(), 2.0)

# outputDiscretization > 1
mesh = mesher3.build(f2d, 2, -1.0, 3.0, 3, True)
print(
    "  outDisc=3: vert=%d simp=%d vol=%.12g"
    % (mesh.getVerticesNumber(), mesh.getSimplicesNumber(), mesh.getVolume())
)
assert mesh.getVerticesNumber() == 16  # 2*2*4
assert mesh.getSimplicesNumber() == 18  # 1*1*3*6
ott.assert_almost_equal(mesh.getVolume(), 2.0)

# 5. 2D non-linear function (original test, kept)
print("=" * 40)
print("2D non-linear function (original test)")
a = [-4.0] * 3
b = [4.0] * 3
f_orig = ot.SymbolicFunction(["x0", "x1"], ["cos(pi_*x0)*sin(pi_*x1)^2"])
f_orig.setName("Paraboloid")
f_orig.setInputDescription([r"$x_0$", r"$x_1$"])
f_orig.setOutputDescription([r"$x_2$"])
inputInterval = ot.Interval([a[0], a[1]], [b[0], b[1]])
inputDiscretization = [100] * 2
mesher = otm.FunctionGraphMesher(inputInterval, inputDiscretization)
print("mesher =", mesher)
outputIndex = 2
outputDiscretization = 1
subGraph = True
mesh = mesher.build(f_orig, outputIndex, a[2], b[2], outputDiscretization, subGraph)
print(mesh)
assert mesh.getDimension() == 3
assert mesh.getName() == "Paraboloid"
assert mesh.getVerticesNumber() == 20402
assert mesh.getSimplicesNumber() == 60000
ott.assert_almost_equal(mesh.getVolume(), 256.0)
assert mesh.getDescription() == [r"$x_0$", r"$x_1$", r"$x_2$"]

# subGraph=False -- supergraph volume equals subgraph volume (symmetric f)
mesh = mesher.build(f_orig, outputIndex, a[2], b[2], outputDiscretization, False)
print("  supergraph volume: %.12g" % mesh.getVolume())
ott.assert_almost_equal(mesh.getVolume(), 256.0)

# outputDiscretization > 1
mesh = mesher.build(f_orig, outputIndex, a[2], b[2], 2, True)
print(
    "  outDisc=2: vert=%d simp=%d vol=%.12g"
    % (mesh.getVerticesNumber(), mesh.getSimplicesNumber(), mesh.getVolume())
)
assert mesh.getVerticesNumber() == 30603  # 101*101*3
assert mesh.getSimplicesNumber() == 120000  # 100*100*2*6
ott.assert_almost_equal(mesh.getVolume(), 256.0)

# outputIndex at 0
mesh = mesher.build(f_orig, 0, a[2], b[2], outputDiscretization, True)
print("  outputIndex=0: vol=%.12g" % mesh.getVolume())
ott.assert_almost_equal(mesh.getVolume(), 256.0)

# 6. Exception tests
print("=" * 40)
print("Exception tests")

with ott.assert_raises(TypeError):
    mesher.build(
        ot.SymbolicFunction(["x0", "x1", "x2"], ["x0"]),
        outputIndex,
        a[2],
        b[2],
        outputDiscretization,
        subGraph,
    )
print("OK: wrong input dim")

with ott.assert_raises(TypeError):
    mesher.build(
        ot.SymbolicFunction(["x0", "x1"], ["x0", "x1"]),
        outputIndex,
        a[2],
        b[2],
        outputDiscretization,
        subGraph,
    )
print("OK: wrong output dim")

with ott.assert_raises(TypeError):
    mesher.build(f_orig, outputIndex, 0.0, 0.0, outputDiscretization, subGraph)
print("OK: equal min/max")

with ott.assert_raises(TypeError):
    mesher.build(f_orig, 5, a[2], b[2], outputDiscretization, subGraph)
print("OK: large outputIndex")

with ott.assert_raises(TypeError):
    mesher.build(f_orig, outputIndex, a[2], b[2], 0, subGraph)
print("OK: zero outputDisc")
