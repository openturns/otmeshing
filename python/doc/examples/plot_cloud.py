"""
Cloud meshing
=============
"""

# %%
import openturns as ot
import openturns.viewer as otv
import otmeshing

# %%
# Using the default triangulation
mesher = otmeshing.CloudMesher()
points = ot.JointDistribution([ot.Uniform()] * 2).getSample(100)
triangulation = mesher.build(points)

# %%
# Plot mesh
graph = triangulation.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Triangulation of sample")
view = otv.View(graph)

# %%
# Now using Delaunay triangulation
mesher = otmeshing.CloudMesher(otmeshing.CloudMesher.DELAUNAY)
triangulation = mesher.build(points)

# %%
# Plot mesh
graph = triangulation.draw()
graph.setLegendPosition("upper left")
graph.setTitle("Delaunay triangulation of sample")
view = otv.View(graph)

# %%
otv.View.ShowAll()
