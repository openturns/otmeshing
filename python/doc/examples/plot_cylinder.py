"""
Cylinder meshing
================
"""

# %%
# In this example we will see how to define generalized cylindric shapes,
# how we to generate their meshes and compute their intersection mesh.
# The Cylinder objects are distinct from their mesh
# and CloudMesher objects will allow to compute their meshes.
# Then knowing their meshes are convex by definition
# we can use the IntersectionMesher to compute the intersection.

# %%
import openturns as ot
import openturns.viewer as otv
import otmeshing as otm

# %%
# Define a 3D cylinder with a disc-shaped base
N = 10
M = 2
dim = 3
print("Build cylinder 1")
xc1 = 0.0
yc1 = 0.0
R1 = 1.0
H1 = 4.0
f1 = ot.SymbolicFunction(
    ["x", "y"], [f"(x-({xc1}))^2+(y-({yc1}))^2"]
)
levelSet1 = ot.LevelSet(f1, ot.LessOrEqual(), R1**2)
base1 = ot.LevelSetMesher([N] * 2).build(
    levelSet1, ot.Interval([xc1 - R1, yc1 - R1], [xc1 + R1, yc1 + R1])
)
extension1 = ot.Interval([-H1 / 2] * (dim - 2), [H1 / 2] * (dim - 2))
injection1 = list(range(2, dim))
cyl1 = otm.Cylinder(base1, extension1, injection1, M)

# %%
# Define an other cylinder in another direction
yc2 = 0.8
zc2 = 0.0
R2 = 0.5
H2 = 4.0
f2 = ot.SymbolicFunction(
    ["y", "z"], [f"(y-({yc2}))^2+(z-({zc2}))^2"]
)
levelSet2 = ot.LevelSet(f2, ot.LessOrEqual(), R2**2)
base2 = ot.LevelSetMesher([N] * 2).build(
    levelSet2, ot.Interval([yc2 - R2, zc2 - R2], [yc2 + R2, zc2 + R2])
)
extension2 = ot.Interval([-H2 / 2] * (dim - 2), [H2 / 2] * (dim - 2))
injection2 = [0] + list(range(3, dim))
cyl2 = otm.Cylinder(base2, extension2, injection2, M)

# %%
# Compute the mesh of both cylinders
method = otm.CloudMesher.BASIC
mesh1 = otm.CloudMesher(method).build(cyl1.getVertices())
mesh2 = otm.CloudMesher(method).build(cyl2.getVertices())

# %%
# Plot cylinder1
shading = False
thetaX, thetaY, thetaZ = 6.1, 2.7, 3.3
drawEdge = True
graph1 = mesh1.draw3D(drawEdge, thetaX, thetaY, thetaZ, shading)
graph1.setLegendPosition("upper left")
graph1.setTitle("Cylinder1 mesh")
view = otv.View(graph1)

# %%
# Plot cylinder2
graph2 = mesh2.draw3D(drawEdge, thetaX, thetaY, thetaZ, shading)
graph2.setLegendPosition("upper left")
graph2.setTitle("Cylinder1 mesh")
view = otv.View(graph2)

# %%
# Compute the intersection, we know both meshes are convex by definition
inter12 = otm.IntersectionMesher().buildConvex([mesh1, mesh2])

# %%
# Plot intersection
graph = inter12.draw3D(drawEdge, thetaX, thetaY, thetaZ, shading)
graph.setLegendPosition("upper left")
graph.setTitle("Cylinder intersection")
view = otv.View(graph)

# %%
otv.View.ShowAll()
