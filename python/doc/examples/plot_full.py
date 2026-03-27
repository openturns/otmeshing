"""
Real-world example
==================
"""

# %%
# This example illustrates the way we can compute the exact intersection of cylinders defining constraints
# in low dimension to get the exact representation of the global admissible set.
# Let's assume that the global parameter space is of dimension :math:`n`, typically :math:`n\in(3,...,10)`.
# Assume that this set is contained inside of a known bounding box
# :math:`[\mathbf{a}, \mathbf{b}]`, and that a point :math:`\mathbf{x}` in :math:`\mathbb{R}^n` is admissible
# if and only if it satisfies a set of constraints :math:`(C_j(\mathbf{x})\leq 0)_{j=1,...,J}`,
# where each of the :math:`C_j` acts only on a low dimensional part :math:`\tilde{\mathbf{x}}\in\mathbb{R}^{n_j}` of :math:`\mathbf{x}`.
# The constraints are defined using one of the following options:
#
#  * The point :math:`\tilde{\mathbf{x}}` must be inside of the convex hull of a set of points in arbitrary dimension;
#  * The point :math:`\tilde{\mathbf{x}}` must be inside of a polygon in dimension 2, given by the collection of its vertices;
#  * The point :math:`\tilde{\mathbf{x}}` must be inside of a given mesh in arbitrary dimension;
#  * The point :math:`\tilde{\mathbf{x}}` must be inside of a sub-graph:
#    there is a function :math:`f_j:\mathbb{R}^{n_j-1}\rightarrow\mathbb{R}`
#    such that for a component :math:`k\in[1,...,n_j], \tilde{\mathbf{x}}_k\leq f(\tilde{\mathbf{x}}_{\sim k})`.
#
# In all these cases, the resulting constraint :math:`C_j` define a cylinder whose **base** is given
# by a mesh in :math:`\mathbb{R}^{n_j}` and whose **extention** is an interval in :math:`\mathbb{R}^{n-n_j}`:
# it is the cartesian product of both objects, with an **injection** :math:`\phi` telling where the components
# of :math:`\tilde{\mathbf{x}}` are positioned into the components of :math:`\mathbf{x}`.

# %%
import math as m
import openturns as ot
import openturns.viewer as otv
import otmeshing as otm
from time import time
ot.RandomGenerator.Reset()

# %%
# We define a set in :math:`\mathbb{R}^3` as the intersection of several constraints acting in :math:`\mathbb{R}^2`
# Here is the global bounding box
a = [-4.0] * 3
b = [4.0] * 3

# %%
# Algorithm used for all the intersection computation
algoInter = otm.IntersectionMesher()

# %%
# Algorithm used for all the convex decomposition computation
# We force the use of simplices decomposition to avoid a bug in the 3D convex
# decomposition (which is correct but awfully slow)
algoDecomp = otm.ConvexDecompositionMesher()
algoDecomp.setUseSimplicesDecomposition(True)

# %%
# Define the first constraint :math:`C_1` as the interior of the convex hull of a set of points in 2D
# :math:`C_1` acts on :math:`(x_0, x_1)`
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
ot.BoundaryMesher().build(C_1.computeMesh()).exportToVTKFile("C_1.vtk")
g1 = ot.BoundaryMesher().build(base1).draw()
cloud = ot.Cloud(points1)
cloud.setPointStyle("fcircle")
g11 = ot.Graph(g1)
g11.add(cloud)
g11
view = otv.View(g11)

# %%
# Define :math:`C_2` as being inside of a polygon given by the vertices of its boundary
# It acts on :math:`(x_0, x_2)`
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
ot.BoundaryMesher().build(C_2.computeMesh()).exportToVTKFile("C_2.vtk")
g2 = ot.BoundaryMesher().build(base2).draw()
cloud = ot.Cloud(points2)
cloud.setPointStyle("fcircle")
g22 = ot.Graph(g2)
g22.add(cloud)
g22
view = otv.View(g22)

# %%
# Define :math:`C_3` as being inside of a mesh built with LevelSetMesher
# It acts on :math:`(x_1, x_2)`
f3 = ot.SymbolicFunction(["x1", "x2"], ["x1^4+x2^3"])
level3 = ot.LevelSet(f3, ot.LessOrEqual(), 4.0)
n3 = 20
base3 = ot.LevelSetMesher([n3] * 2).build(level3, ot.Interval([a[1], a[2]], [b[1], b[2]]))
base3.setName(r"$C_3$")
base3.setDescription([r"$x_1$", r"$x_2$"])
extension3 = ot.Interval(a[0], b[0])
injection3 = [0]
N3 = 1
C_3 = otm.Cylinder(base3, extension3, injection3, N3)
ot.BoundaryMesher().build(C_3.computeMesh()).exportToVTKFile("C_3.vtk")
g3 = ot.BoundaryMesher().build(base3).draw()
g3
view = otv.View(g3)

# %%
# Now the admissible domain
meshAllCylinders = algoInter.buildCylinder([C_1, C_2, C_3])
ot.BoundaryMesher().build(meshAllCylinders).exportToVTKFile("cylindersIntersection.vtk")
meshAllCylindersConvexParts = algoDecomp.build(meshAllCylinders)
print("Number of convex parts=", len(meshAllCylindersConvexParts))

# %%
# Create a uniform distribution over mesh and sample it
distribution = ot.UniformOverMesh(meshAllCylinders)
size = 1000
sample = distribution.getSample(size)

# %%
# Check if all the constraints are satisfied
grid = ot.GridLayout(1, 3)
g1c = ot.Graph(g1)
g1c.add(ot.Cloud(sample.getMarginal([0, 1])))
grid.setGraph(0, 0, g1c)
g2c = ot.Graph(g2)
g2c.add(ot.Cloud(sample.getMarginal([0, 2])))
grid.setGraph(0, 1, g2c)
g3c = ot.Graph(g3)
g3c.add(ot.Cloud(sample.getMarginal([1, 2])))
grid.setGraph(0, 2, g3c)
grid
view = otv.View(grid)

# %%
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
ot.BoundaryMesher().build(mesh).exportToVTKFile("func_graph.vtk")

meshConvexParts = algoDecomp.build(mesh)
print("Number of convex parts=", len(meshConvexParts))

# %%
t0 = time()
globalMesh = algoInter.build([mesh, meshAllCylinders])
t1 = time()
ot.BoundaryMesher().build(globalMesh).exportToVTKFile("global.vtk")
print("t=", t1 - t0, "s")

# %%
# Create a uniform distribution over mesh and sample it
distribution = ot.UniformOverMesh(globalMesh)
size = 100000
sample = distribution.getSample(size)
ot.Mesh(sample).exportToVTKFile("Global_sample.vtk")
# Extract only the 1000 first points to check things in 2D
sample = sample[:1000]
# Check if all the constraints are satisfied
grid = ot.GridLayout(1, 3)
g1g = ot.Graph(g1)
g1g.add(ot.Cloud(sample.getMarginal([0, 1])))
grid.setGraph(0, 0, g1g)
g2g = ot.Graph(g2)
g2g.add(ot.Cloud(sample.getMarginal([0, 2])))
grid.setGraph(0, 1, g2g)
g3g = ot.Graph(g3)
g3g.add(ot.Cloud(sample.getMarginal([1, 2])))
grid.setGraph(0, 2, g3g)
grid
view = otv.View(grid)

# %%
otv.View.ShowAll()
