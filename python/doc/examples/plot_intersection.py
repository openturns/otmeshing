"""
Intersection meshing
====================
"""

# %%
# In this example we will see how to define a mesh from the intersection of several meshes.

# %%
import math
import openturns as ot
import openturns.viewer as otv
import otmeshing as otm

# %%
# First mesh: a box
mesh1 = ot.IntervalMesher([1] * 2).build(ot.Interval([-1.2] * 2, [1.2] * 2))

# %%
# Plot first mesh
graph = mesh1.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Mesh 1")
view = otv.View(graph)

# %%
# Second mesh: an arch
mesh2 = ot.IntervalMesher([10, 10]).build(ot.Interval([1.0, 0.0], [2.0, math.pi]))
f = ot.SymbolicFunction(["r", "theta"], ["r*cos(theta)", "r*sin(theta)"])
mesh2.setVertices(f(mesh2.getVertices()))

# %%
# Plot second mesh
graph = mesh2.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Mesh 2")
view = otv.View(graph)

# %%
# Compute the intersection
mesher = otm.IntersectionMesher()
intersection = mesher.build([mesh1, mesh2])

# %%
# Plot union
graph = intersection.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Mesh intersection")
view = otv.View(graph)

# %%
otv.View.ShowAll()
