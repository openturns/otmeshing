"""
Union meshing
=============
"""

# %%
# In this example we will see how to define a mesh from the union of several meshes.

# %%
import openturns as ot
import openturns.viewer as otv
import otmeshing as otm

# %%
# Define two non-overlapping meshes
dim = 2
mesh1 = ot.IntervalMesher([1] * dim).build(ot.Interval(dim))
mesh2 = ot.IntervalMesher([1] * dim).build(ot.Interval([2.0] * dim, [3.0] * dim))

# %%
# Compute the mesh of both cylinders
mesher = otm.UnionMesher()
union = mesher.build([mesh1, mesh2])

# %%
# Plot union
graph = union.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Mesh union")
view = otv.View(graph)

# %%
otv.View.ShowAll()
