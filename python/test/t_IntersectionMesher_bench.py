import math as m
import openturns as ot
import otmeshing as otm
from time import time

algoInter = otm.IntersectionMesher()
algoDecomp = otm.ConvexDecompositionMesher()

# the default is simplicial decomposition but with coacd enabled
# (instead of cgal) the decomposition step becomes interesting
if otm.ConvexDecompositionMesher.HasFeature("coacd"):
    algoDecomp.setUseSimplicesDecomposition(False)
    algoInter.setUseSimplicesDecomposition(False)

# C1
a = [-4.0] * 3
b = [4.0] * 3
points1 = ot.Normal(2).getSample(20)
# use ConvexHullMesher to exclude internal points
hull1 = otm.ConvexHullMesher().build(points1).getVertices()
base1 = otm.CloudMesher().build(hull1)
base1.setName(r"$C_1$")
base1.setDescription([r"$x_0$", r"$x_1$"])
extension1 = ot.Interval(a[2], b[2])
injection1 = [2]
N1 = 1
C_1 = otm.Cylinder(base1, extension1, injection1, N1)

# C2
t = 0.0
dt = 1e-6
eps = 2e-1
points2 = [ot.Point([1, 2])]
# [cos(t), sin(t)+2cos(2t), t=0..2pi]
while t + dt < 2 * m.pi:
    p = ot.Point([m.cos(t + dt), m.sin(t + dt) + 2 * m.cos(2 * (t + dt))])
    dp = (p - points2[-1]).norm()
    while dp < 0.5 * eps:
        dt *= 1.1
        p = ot.Point([m.cos(t + dt), m.sin(t + dt) + 2 * m.cos(2 * (t + dt))])
        dp = (p - points2[-1]).norm()
    while dp > 2.0 * eps:
        dt /= 1.1
        p = ot.Point([m.cos(t + dt), m.sin(t + dt) + 2 * m.cos(2 * (t + dt))])
        dp = (p - points2[-1]).norm()
    points2.append(p)
    t += dt
base2 = otm.PolygonMesher().build(points2)
base2.setName(r"$C_2$")
base2.setDescription([r"$x_0$", r"$x_2$"])
extension2 = ot.Interval(a[1], b[1])
injection2 = [1]
N2 = 1
C_2 = otm.Cylinder(base2, extension2, injection2, N2)


# C3
f3 = ot.SymbolicFunction(["x1", "x2"], ["x1^4+x2^3"])
level3 = ot.LevelSet(f3, ot.LessOrEqual(), 4.0)
n3 = 8
base3 = ot.LevelSetMesher([n3] * 2).build(level3, ot.Interval([a[1], a[2]], [b[1], b[2]]))
base3.setName(r"$C_3$")
base3.setDescription([r"$x_1$", r"$x_2$"])
extension3 = ot.Interval(a[0], b[0])
injection3 = [0]
N3 = 1
C_3 = otm.Cylinder(base3, extension3, injection3, N3)

for c in [C_1, C_2, C_3]:
    print(f"convex ? {c.isConvex()}")

# domain=cylinders intersection, keep them as list of convexes for later
t0 = time()
convexPiecesAllCylinders = [algoInter.buildConvexSample([c.getVertices() for c in [C_1, C_2, C_3]])]
t1 = time()
print(f"buildCylinderConvex t={t1 - t0} s")
print(f"convexPiecesAllCylinders size={len(convexPiecesAllCylinders)}")

# FGM
f = ot.SymbolicFunction(["x0", "x1"], ["1+2*cos(pi_*x0/2)*sin(pi_*x1/2)^2"])
f.setName("Paraboloid")
f.setInputDescription([r"$x_0$", r"$x_1$"])
f.setOutputDescription([r"$x_2$"])
inputInterval = ot.Interval([a[0], a[1]], [b[0], b[1]])
inputDiscretization = [41] * 2
outputDimension = 2
outputDiscretization = 1
mesher = otm.FunctionGraphMesher(inputInterval, inputDiscretization)
mesh = mesher.build(f, outputDimension, a[2], b[2], 1)

# decompose FGM
t0 = time()
meshConvexParts = algoDecomp.build(mesh)
t1 = time()
print("Number of convex parts=", len(meshConvexParts))
print(f"meshConvexParts t={t1 - t0} s")

# intersection cylinders / FGM
# pass directly the cylinder intersection as decomposition of convexes to avoid
# the final assembly step yielding incompatible topology on shared faces
# because of the independent triangulations of each convex component
# (every edge must be shared by exactly two triangles with opposite orientation)
# else CoACD would throw "The mesh is not a 2-manifold!".
t0 = time()
globalMesh = algoInter.buildConvex([mesh] + [otm.CloudMesher().build(c) for c in convexPiecesAllCylinders])
t1 = time()
print(f"build t={t1 - t0} s")
